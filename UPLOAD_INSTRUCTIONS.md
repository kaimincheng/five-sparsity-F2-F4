# Upload checklist (for the authors)

1. Create or open the repository `kaimincheng/five-sparsity-F2-F4`.
2. **Extract** the prepared ZIP on your computer. Do not upload the ZIP as
   the only repository file: GitHub does not automatically unpack it.
3. On the repository's **Code** tab, choose **Add file → Upload files**.
   Upload all extracted files and the `data/` and `logs/` directories,
   preserving their relative paths. If a README already exists, use this
   `README.md` as the replacement after reviewing its contents.
4. Commit your changes, e.g. `Add exact verification certificates C1-C9`.
5. Confirm that all source files, the CSV, and the log are visible on GitHub.
6. On a local machine with Python 3 and a C++17 compiler, run
   `bash run_all.sh`; it should end with `ALL CHECKS C1--C9 PASSED`.
7. Compare the paper's Section 7 with the uploaded files and results.
8. For a submission-ready snapshot, create a release (e.g. `v1.0.0`)
   or record the exact commit hash, and put that stable reference in the paper.

If you use GitHub in your browser, you must drag the **files after extracting**,
not just the ZIP. If the browser does not keep folder hierarchy, upload files
in `data/` and `logs/` using those paths via the GitHub file editor or Git.

GitHub repository: https://github.com/kaimincheng/five-sparsity-F2-F4
