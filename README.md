## How to get and run the code

This page is hosted on GitHub. You do **not** need to know GitHub or Git to use the code.

A quick note about GitHub buttons:

- **Copy** copies text to your clipboard. It does not download a runnable file.
- **Raw** shows the plain file contents in the browser.
- **Download raw file** downloads the actual file to your computer.
- **Download ZIP** downloads the whole repository.

At the moment, the main runnable file in this repository is:

| File | What it does | How to run it |
|------|--------------|---------------|
| [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb) | Computes the third-order cooperation index χ(i,j,k) for two uploaded PDB files | Run in Google Colab |

---

## Option 1 — Run the notebook in Google Colab

This is the easiest option for readers who do not use GitHub.

1. Open `compute_chi_ijk.ipynb` on GitHub.
2. Click **Raw** or **Download raw file** if you want to save the notebook.
3. Open [Google Colab](https://colab.research.google.com/).
4. In Colab, choose **File → Upload notebook**.
5. Upload `compute_chi_ijk.ipynb`.
6. Run the cells from top to bottom.
7. When prompted, upload exactly two `.pdb` files.

The notebook will compute χ(i,j,k) and create output CSV files and a ZIP file for download.

---

## Option 2 — Download the whole repository

Use this if you want a local copy of everything.

1. Go to the repository main page.
2. Click the green **`<> Code`** button.
3. Choose **Download ZIP**.
4. Unzip the downloaded file.
5. Open the notebook file `compute_chi_ijk.ipynb`.

You can then upload the notebook to Google Colab as described above.

---

## Option 3 — Use Git, if you already know Git

If you already use Git, you can clone the repository:

```bash
git clone https://github.com/fatmasenguler/laplacian-minor-hierarchy.git
cd laplacian-minor-hierarchy
