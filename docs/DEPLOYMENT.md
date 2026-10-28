# Deployment guide

How to put the Streamlit interface online, and how to check it before you do.

> **Status:** deployed. The live app is
> <https://dnasequencetoolkit.streamlit.app/>, built from the `main` branch with
> `app.py` as the entry point. Pushing to `main` triggers a redeploy.

## 1. Check it runs locally first

```bash
pip install -r requirements.txt
streamlit run app.py
```

Streamlit prints a local URL (usually <http://localhost:8501>). Open it and
confirm:

- the sample sequence loads and the seven sections populate;
- pasting `ATGCX` produces an error message rather than a crash;
- uploading a two-record FASTA gives you a record selector;
- a CSV download actually downloads.

If any of those fail locally, they will fail in the cloud too.

## 2. Push to GitHub

The app must live in a **public** repository for the free tier.

```bash
git add .
git commit -m "Add Streamlit interface and documentation"
git push origin main
```

## 3. Deploy on Streamlit Community Cloud

1. Go to <https://share.streamlit.io> and sign in with GitHub.
2. Click **New app**.
3. Choose **Deploy a public app from GitHub**.
4. Fill in:
   - **Repository:** `manojkumarmansingh/dna-sequence-toolkit`
   - **Branch:** `main`
   - **Main file path:** `app.py`
5. Optionally set a custom subdomain.
6. Click **Deploy**.

Streamlit installs the packages listed in `requirements.txt` and runs `app.py`.
First build usually takes a couple of minutes.

## 4. How dependencies resolve

| File | Used by | Contents |
| --- | --- | --- |
| `requirements.txt` | Streamlit Cloud | `streamlit`, `pandas` |
| `requirements-dev.txt` | local testing / CI | `pytest` |
| `pyproject.toml` | `pip install -e .` | the `dna_toolkit` package itself |

`dna_toolkit` has no third-party runtime dependencies; `requirements.txt`
exists only for the web interface.

## 5. Common problems

- **`ModuleNotFoundError: dna_toolkit`**: the app imports the package from the
  repository root. Deploy from the repository root, not a subfolder, and make
  sure `app.py` sits at the top level.
- **`ModuleNotFoundError: examples`**: `app.py` already falls back gracefully,
  but the sample picker needs `examples/samples.py` in the repository.
- **App sleeps**: free-tier apps sleep after a period of inactivity; the next
  visitor triggers a cold start.

## 6. Screenshots

Capture these once the app is running, and save them to `docs/screenshots/`:

| File | What to capture |
| --- | --- |
| `overview.png` | Overview tab with the sample sequence loaded |
| `orfs.png` | ORF table with several rows |
| `restriction.png` | Restriction-site table |
| `error.png` | The validation message after entering an invalid character |

Then add to the README:

```markdown
![Overview](docs/screenshots/overview.png)
```

Keep images under ~500 KB each; PNG or WebP both work on GitHub.

## 7. Release notes

When you are happy with the deployment, publish a GitHub release:

1. **Releases** → **Draft a new release**.
2. Tag: `v0.2.0` (create it on publish).
3. Title: `v0.2.0: Streamlit interface and hardened core`.
4. Paste the `[0.2.0]` section from [CHANGELOG.md](../CHANGELOG.md).

Suggested one-paragraph summary:

> Adds a Streamlit web interface with CSV/FASTA/JSON downloads, configurable ORF
> start codons, circular restriction mapping, record-level FASTA validation and
> a 36-test edge-case suite. The public API is unchanged; all new behaviour is
> opt-in.
