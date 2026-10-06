# De Bruijn Graph Genome Construction

A web-based application for genome construction using the De Bruijn Graph approach. The application accepts DNA sequences, generates k-mers, constructs a directed De Bruijn graph, calculates frequency-based edge weights, analyzes the graph, finds a suitable path, and reconstructs the DNA sequence.

## Features

- Selection of 10 provided DNA sequences: S1–S10
- Custom DNA sequence input
- TXT/FASTA file upload
- DNA sequence validation
- User-defined k-mer size
- Automatic k-mer generation
- De Bruijn graph construction
- `(k-1)`-mers represented as nodes
- k-mers represented as directed edges
- Frequency-based edge weights
- Edge weight normalization
- In-degree and out-degree calculation
- Eulerian path/circuit analysis
- Eulerian path using Hierholzer's algorithm
- Frequency-weighted best-path analysis
- DNA sequence reconstruction
- Interactive graph visualization
- k-mer and weight tables
- Graph statistics
- All calculations based on the current input

## Technologies Used

- Python
- FastAPI
- Uvicorn
- HTML5
- CSS3
- JavaScript
- SVG

## De Bruijn Graph Concept

A De Bruijn graph represents DNA sequences using overlapping k-mers.

For example, consider the DNA sequence `ATGATGATGCAT` with `k = 3`.

The generated k-mers are:

`ATG`, `TGA`, `GAT`, `ATG`, `TGA`, `GAT`, `ATG`, `TGC`, `GCA`, `CAT`

For a De Bruijn graph:

- Nodes represent `(k-1)`-mers.
- Edges represent k-mers.
- The prefix of a k-mer becomes the source node.
- The suffix of a k-mer becomes the destination node.

For example, the k-mer `ATG` produces the edge:

`AT → TG`

## Mathematical Model

Let the input DNA sequence be `S` with length `n`.

For a selected k-mer size `k`, the number of k-mers is:

`m = n - k + 1`

Each k-mer is represented as:

`Kᵢ = sᵢsᵢ₊₁...sᵢ₊ₖ₋₁`

For every k-mer:

`prefix(K) = first (k-1) characters`

`suffix(K) = last (k-1) characters`

The De Bruijn graph is represented as:

`G = (V, E)`

where:

- `V` is the set of unique `(k-1)`-mers.
- `E` is the set of k-mers.

Each edge is constructed as:

`prefix(K) → suffix(K)`

## Edge Weight Calculation

Each edge is assigned a weight according to the frequency of its corresponding k-mer in the current input sequence.

The formula is:

`W(eK) = f(K)`

where `f(K)` represents the number of occurrences of k-mer `K`.

For example, for the sequence `ATGATGATGCAT`:

| K-mer | Frequency | Edge Weight |
|---|---:|---:|
| ATG | 3 | 3 |
| TGA | 2 | 2 |
| GAT | 2 | 2 |
| TGC | 1 | 1 |
| GCA | 1 | 1 |
| CAT | 1 | 1 |

The weights are recalculated whenever the input sequence changes.

## Weight Normalization

The maximum edge weight is calculated as:

`Wmax = max(W(e))`

The normalized edge weight is:

`Wnorm(e) = W(e) / Wmax`

For example, when `Wmax = 3`:

| K-mer | Weight | Normalized Weight |
|---|---:|---:|
| ATG | 3 | 1.000 |
| TGA | 2 | 0.667 |
| GAT | 2 | 0.667 |
| TGC | 1 | 0.333 |
| GCA | 1 | 0.333 |
| CAT | 1 | 0.333 |

## Graph Analysis

After constructing the graph, the application calculates:

- Number of k-mers
- Number of unique k-mers
- Number of nodes
- Number of edges
- k-mer frequencies
- Edge weights
- Normalized edge weights
- Node in-degree
- Node out-degree
- Eulerian path/circuit conditions

## Path Finding

### Eulerian Path

When the graph satisfies the required Eulerian conditions, Hierholzer's algorithm is used.

The objective is to traverse every edge exactly once.

### Weighted Best Path

The application also performs weighted path analysis using the calculated edge weights.

The score of a path `P` is:

`Score(P) = Σ W(e)`

A path containing higher-frequency k-mers receives a higher score.

The weighted approach is particularly useful when multiple possible routes exist in the graph.

## Genome Reconstruction

After obtaining an ordered sequence of k-mers, the DNA sequence is reconstructed.

The first k-mer is taken completely. For every subsequent k-mer, only its final nucleotide is appended.

For example:

`ATG → TGA → GAT → ATG`

produces:

`ATGATGATG`

## Dataset

The application contains 10 sequencing inputs:

- S1
- S2
- S3
- S4
- S5
- S6
- S7
- S8
- S9
- S10

The user can select any of these sequences from the dropdown.

The application also supports:

- Custom DNA sequence input
- TXT file input
- FASTA file input

All calculations are performed using the input currently selected or provided by the user.

## Application Workflow

Input DNA Sequence → Validate Sequence → Select k → Generate k-mers → Calculate k-mer Frequency → Calculate Edge Weights → Normalize Weights → Construct De Bruijn Graph → Analyze Graph → Find Path → Reconstruct DNA Sequence → Display Results

## Project Structure

The project can be organized as follows:

- `app.py` — FastAPI backend and application logic
- `sequences.txt` — DNA sequence dataset
- `requirements.txt` — Python dependencies
- `README.md` — Project documentation
- `.gitignore` — Git ignored files
- `templates/index.html` — Web interface
- `static/style.css` — User interface styling
- `static/script.js` — Frontend functionality

## Installation

### Clone the Repository

`git clone https://github.com/YOUR_USERNAME/DeBruijn-Genome-Construction.git`

### Open the Project Directory

`cd DeBruijn-Genome-Construction`

### Create a Virtual Environment

`python -m venv venv`

### Activate the Virtual Environment

For Windows:

`venv\Scripts\activate`

For Linux/macOS:

`source venv/bin/activate`

### Install Dependencies

`pip install -r requirements.txt`

## Running the Application

Run the application using:

`python app.py`

Alternatively, use Uvicorn:

`python -m uvicorn app:app --reload`

After starting the server, open:

`http://127.0.0.1:8000`

in a web browser.

## Example

Consider the input sequence:

`ATGATGATGCAT`

with:

`k = 3`

The application generates 10 k-mers.

The graph contains nodes such as:

`AT`, `TG`, `GA`, `GC`, `CA`

and edges such as:

`ATG`, `TGA`, `GAT`, `TGC`, `GCA`, `CAT`

The edge weights are calculated from their frequencies in the input sequence.

The resulting graph and path information are then displayed in the web interface.

## Output

The application displays:

- Selected sequence
- Input DNA sequence
- k value
- Generated k-mers
- k-mer frequencies
- Edge weights
- Normalized weights
- Graph nodes
- Graph edges
- Node degrees
- Eulerian status
- Selected/best path
- Path score
- Reconstructed DNA sequence
- Graph visualization
- Graph statistics

## Future Scope

The project can be extended with:

- Multiple-read genome assembly
- Contig generation
- Read error correction
- Coverage-based weighting
- Quality-score-based weighting
- Paired-end read support
- Larger genome datasets
- Improved graph-layout algorithms
- Assembly quality metrics
- FASTA export
- Advanced genome assembly techniques

## Author

**Anish Dey**

M.Tech in Artificial Intelligence and Machine Learning with Pedagogy

Department of Computer Science and Engineering

Session: 2026–27

## License

This project is developed for academic and educational purposes.