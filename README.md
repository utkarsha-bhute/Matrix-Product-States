# Matrix Product States (MPS) for Qudit Systems

This repository contains Python code to construct various canonical forms of Matrix Product States (MPS) for **qudit systems** (i.e., local dimension *d* > 2). The implementations include:

- Left-canonical MPS (LCMPS)
- Right-canonical MPS (RCMPS)
- Mixed-canonical MPS (MCMPS)
- Functions for overlap, expectation values, and reduced density matrices (RDM)

---

## 📁 Repository Structure

- `mps_qudit.py` – Core implementation of MPS construction and related functions
- `test_file.txt` – Output file where matrices are saved
- `README.md` – Project documentation (this file)

---

## ⚙️ Features

- Normalize input quantum state vectors
- Construct left, right, and mixed-canonical MPS
- Compute:
  - Overlap ⟨ψ₁|ψ₂⟩ between two MPS states
  - Expectation value ⟨ψ|O|ψ⟩ of local operators
  - Matrix element ⟨ψ₂|O|ψ₁⟩
  - Reduced density matrices (RDM) for bipartitions

---

## 🧠 Background

Matrix Product States (MPS) are tensor network representations useful for simulating quantum many-body systems, especially those with limited entanglement. This code explores **qudit** systems with customizable local dimensions and system sizes.

---

## 🐍 Requirements

- Python 3.x
- NumPy

Install NumPy using:

```bash
pip install numpy
```

---

## 🚀 Usage
## 1. Left-Canonical MPS
```bash
from mps_qudit import normalize, leftCanonicalMPS

L = 3        # Number of sites
d = 3        # Local dimension (qudit)
psir = [0]*d**L
psir[0], psir[13], psir[26] = 1, 1, 1

vin = normalize(psir)
matricesA = leftCanonicalMPS(d, L, vin.reshape((1, -1)), 1)
```

## 2. Right-Canonical MPS
```bash
from mps_qudit import normalize, rightCanonicalMPS

L = 4
d = 2
psir = [0]*d**L
psir[0], psir[15] = 1, 1

vin = normalize(psir)
matricesB = rightCanonicalMPS(d, L, vin.reshape((1, -1)), 1)
```

## 3. Mixed-Canonical MPS
```bash
from mps_qudit import normalize, mixedCanonicalMPS

L = 3
d = 3
p = 2  # Canonical center
psir = [0]*d**L
psir[0], psir[13], psir[26] = 1, 1, 1

vin = normalize(psir)
matrices = mixedCanonicalMPS(d, L, vin.reshape((1, -1)), 1, p)
```
---

## 📄 Output Format

The MPS matrices are written to ``` test_file.txt ``` in the following format:
```bash
A Matrices:
A0 =
[[...]]
A1 =
[[...]]

B Matrices:
B0 =
[[...]]
B1 =
[[...]]
```

---

## 👤 Author
Utkarsha Bhute

📫 Feel free to reach out for feedback, questions, or collaboration.

