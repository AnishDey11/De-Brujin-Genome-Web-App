from collections import Counter, defaultdict, deque
from pathlib import Path
import re

from fastapi import FastAPI, File, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "sequences.txt"

app = FastAPI(title="Weighted De Bruijn Graph - Genome Construction")
app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")


def clean_dna(text: str) -> str:
    return re.sub(r"[^ACGT]", "", str(text).upper())


def parse_sequences(text: str) -> dict[str, str]:
    """Parse normal multi-record FASTA or the supplied one-header/many-read TXT format."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return {}

    header_count = sum(line.startswith(">") for line in lines)

    if header_count > 1:
        result = {}
        current = None
        for line in lines:
            if line.startswith(">"):
                raw = line[1:].strip()
                base = raw or f"S{len(result) + 1}"
                name = base
                n = 2
                while name in result:
                    name = f"{base}_{n}"
                    n += 1
                result[name] = ""
                current = name
            elif current:
                result[current] += clean_dna(line)
        return {k: v for k, v in result.items() if v}

    # The supplied file has one header followed by ten independent read lines.
    dna_lines = [clean_dna(line) for line in lines if not line.startswith(">")]
    dna_lines = [seq for seq in dna_lines if seq]
    return {f"S{i + 1}": seq for i, seq in enumerate(dna_lines)}


def load_sequences_from_disk() -> dict[str, str]:
    if DATA_FILE.exists():
        parsed = parse_sequences(DATA_FILE.read_text(encoding="utf-8"))
        if parsed:
            return parsed
    return {}


SEQUENCES = load_sequences_from_disk()


def kmer_counts(sequence: str, k: int) -> Counter:
    sequence = clean_dna(sequence)
    if k < 1 or len(sequence) < k:
        return Counter()
    return Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def build_graph(sequence: str, k: int, counts: Counter) -> dict:
    sequence = clean_dna(sequence)
    edges = []
    adjacency = defaultdict(list)
    indeg = Counter()
    outdeg = Counter()

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        u = kmer[:-1]
        v = kmer[1:]
        edge = {
            "id": i,
            "kmer": kmer,
            "from": u,
            "to": v,
            "weight": counts[kmer],
        }
        edges.append(edge)
        adjacency[u].append(i)
        outdeg[u] += 1
        indeg[v] += 1

    nodes = sorted(set(indeg) | set(outdeg))
    return {
        "nodes": nodes,
        "edges": edges,
        "adj": {k: list(v) for k, v in adjacency.items()},
        "indeg": dict(indeg),
        "outdeg": dict(outdeg),
    }


def degree_status(graph: dict) -> tuple[bool, list[str], list[str]]:
    plus = [n for n in graph["nodes"] if graph["outdeg"].get(n, 0) - graph["indeg"].get(n, 0) == 1]
    minus = [n for n in graph["nodes"] if graph["indeg"].get(n, 0) - graph["outdeg"].get(n, 0) == 1]
    bad = [n for n in graph["nodes"] if abs(graph["outdeg"].get(n, 0) - graph["indeg"].get(n, 0)) > 1]
    valid = not bad and ((len(plus) == 1 and len(minus) == 1) or (len(plus) == 0 and len(minus) == 0))
    return valid, plus, minus


def eulerian_path(graph: dict) -> dict:
    edges = graph["edges"]
    if not edges:
        return {"is_eulerian": False, "edge_ids": [], "nodes": [], "score": 0}

    valid_degree, plus, minus = degree_status(graph)
    if plus:
        start = plus[0]
    else:
        start = next((n for n in graph["nodes"] if graph["outdeg"].get(n, 0) > 0), graph["nodes"][0])

    # Hierholzer with edge IDs retained directly. This fixes duplicate-edge bookkeeping.
    local = {node: list(ids) for node, ids in graph["adj"].items()}
    stack_nodes = [start]
    stack_edges = []
    reverse_nodes = []
    reverse_edges = []

    while stack_nodes:
        u = stack_nodes[-1]
        if local.get(u):
            eid = local[u].pop()
            stack_nodes.append(edges[eid]["to"])
            stack_edges.append(eid)
        else:
            reverse_nodes.append(stack_nodes.pop())
            if stack_edges:
                reverse_edges.append(stack_edges.pop())

    path_nodes = list(reversed(reverse_nodes))
    path_edges = list(reversed(reverse_edges))
    used_all = len(path_edges) == len(edges)
    is_eulerian = valid_degree and used_all and len(path_nodes) == len(edges) + 1
    score = sum(edges[eid]["weight"] for eid in path_edges)
    return {"is_eulerian": is_eulerian, "edge_ids": path_edges, "nodes": path_nodes, "score": score}


def weighted_best_path(graph: dict) -> dict:
    """Greedy weighted trail with local backtracking fallback.

    Every edge has a positive frequency. If an Eulerian traversal exists, the
    highest possible score among full-edge traversals is the sum of all edge
    weights, so the Eulerian path is used as the weighted best path. Otherwise,
    a bounded best-score trail is found by exploring high-weight unused edges.
    """
    edges = graph["edges"]
    if not edges:
        return {"edge_ids": [], "nodes": [], "score": 0, "method": "No edges"}

    euler = eulerian_path(graph)
    if euler["is_eulerian"]:
        return {
            "edge_ids": euler["edge_ids"],
            "nodes": euler["nodes"],
            "score": euler["score"],
            "method": "Eulerian traversal with frequency weights",
        }

    # Bounded DFS/beam for non-Eulerian input. Keep the top states at each depth.
    beam_width = 80
    max_depth = min(len(edges), 250)
    states = []
    for node in graph["nodes"]:
        if graph["adj"].get(node):
            states.append((node, 0, [], [node], frozenset()))

    best = max(states, key=lambda x: x[1], default=None)
    for _ in range(max_depth):
        next_states = []
        for node, score, eids, nodes, used in states:
            candidates = sorted(graph["adj"].get(node, []), key=lambda eid: edges[eid]["weight"], reverse=True)
            for eid in candidates:
                if eid in used:
                    continue
                edge = edges[eid]
                nscore = score + edge["weight"]
                nused = used | {eid}
                next_states.append((edge["to"], nscore, eids + [eid], nodes + [edge["to"]], nused))

        if not next_states:
            break
        next_states.sort(key=lambda x: (x[1], len(x[2])), reverse=True)
        states = next_states[:beam_width]
        candidate = states[0]
        if best is None or (candidate[1], len(candidate[2])) > (best[1], len(best[2])):
            best = candidate

    if best is None:
        return {"edge_ids": [], "nodes": [], "score": 0, "method": "No valid trail"}
    return {"edge_ids": best[2], "nodes": best[3], "score": best[1], "method": "Frequency-weighted best trail (bounded search)"}


def reconstruct(edges: list[dict], edge_ids: list[int]) -> str:
    if not edge_ids:
        return ""
    result = edges[edge_ids[0]]["kmer"]
    for eid in edge_ids[1:]:
        result += edges[eid]["kmer"][-1]
    return result


def browser_graph(graph: dict, max_edges: int = 250):
    edges = graph["edges"][:max_edges]
    used_nodes = set()
    for edge in edges:
        used_nodes.add(edge["from"])
        used_nodes.add(edge["to"])
    max_weight = max((e["weight"] for e in graph["edges"]), default=1)
    nodes = [
        {"id": n, "in": graph["indeg"].get(n, 0), "out": graph["outdeg"].get(n, 0)}
        for n in used_nodes
    ]
    browser_edges = [{**e, "normalized_weight": round(e["weight"] / max_weight, 4)} for e in edges]
    return nodes, browser_edges


def analyze_sequence(sequence: str, sequence_id: str, k: int, source_label: str) -> dict:
    sequence = clean_dna(sequence)
    if not sequence:
        raise ValueError("No valid DNA bases were found. Use only A, C, G and T.")
    if k < 2:
        raise ValueError("k must be at least 2.")
    if len(sequence) < k:
        raise ValueError(f"The input contains {len(sequence)} bases, but k={k}. Choose a smaller k.")

    # IMPORTANT: all calculations use ONLY the input sequence being analyzed.
    counts = kmer_counts(sequence, k)
    graph = build_graph(sequence, k, counts)
    euler = eulerian_path(graph)
    best = weighted_best_path(graph)
    nodes, browser_edges = browser_graph(graph)
    max_weight = max((e["weight"] for e in graph["edges"]), default=0)

    all_edges = [
        {**e, "normalized_weight": round(e["weight"] / max(max_weight, 1), 4)}
        for e in graph["edges"]
    ]

    common = [{"kmer": kmer, "frequency": freq} for kmer, freq in counts.most_common(100)]
    return {
        "sequence_id": sequence_id,
        "source": source_label,
        "sequence": sequence,
        "sequence_length": len(sequence),
        "k": k,
        "nodes_count": len(graph["nodes"]),
        "edges_count": len(graph["edges"]),
        "unique_kmers": len(counts),
        "max_weight": max_weight,
        "eulerian": euler["is_eulerian"],
        "euler_score": euler["score"],
        "euler_edges_used": len(euler["edge_ids"]),
        "euler_reconstruction": reconstruct(graph["edges"], euler["edge_ids"]),
        "best_score": best["score"],
        "best_edges_used": len(best["edge_ids"]),
        "best_path_nodes": best["nodes"],
        "best_reconstruction": reconstruct(graph["edges"], best["edge_ids"]),
        "best_method": best["method"],
        "nodes": nodes,
        "edges": browser_edges,
        "all_edges": all_edges,
        "common_kmers": common,
    }


@app.get("/", response_class=HTMLResponse)
async def home():
    return HTMLResponse((BASE_DIR / "templates" / "index.html").read_text(encoding="utf-8"))


@app.get("/api/sequences")
async def get_sequences():
    return {"sequences": {name: len(seq) for name, seq in SEQUENCES.items()}}


@app.get("/api/sequence-text")
async def get_sequence_text(sequence_id: str = "S1"):
    if sequence_id not in SEQUENCES:
        return JSONResponse({"error": "Sequence not found."}, status_code=404)
    return {"sequence_id": sequence_id, "sequence": SEQUENCES[sequence_id]}


@app.get("/api/analyze")
async def analyze(sequence_id: str = "S1", k: int = 21):
    if sequence_id not in SEQUENCES:
        return JSONResponse({"error": "Sequence not found."}, status_code=404)
    try:
        return analyze_sequence(SEQUENCES[sequence_id], sequence_id, k, f"Dataset read {sequence_id}")
    except ValueError as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)


@app.post("/api/analyze-custom")
async def analyze_custom(payload: dict):
    try:
        sequence = str(payload.get("sequence", ""))
        k = int(payload.get("k", 21))
        return analyze_sequence(sequence, "CUSTOM", k, "Typed DNA sequence")
    except (ValueError, TypeError) as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)


@app.post("/api/upload")
async def upload(file: UploadFile = File(...)):
    global SEQUENCES
    content = await file.read()
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        text = content.decode("latin-1")

    parsed = parse_sequences(text)
    if not parsed:
        return JSONResponse({"error": "No valid DNA sequences were found in the file."}, status_code=400)

    SEQUENCES = parsed
    DATA_FILE.write_text(text, encoding="utf-8")
    return {"success": True, "sequences": {name: len(seq) for name, seq in SEQUENCES.items()}}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app:app", host="127.0.0.1", port=8000, reload=True)
