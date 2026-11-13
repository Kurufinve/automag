# Prophet model checkpoints

The two checkpoints of [Prophet-V1-Mag](https://github.com/kairosmaterial/prophet)
that `calculator: prophet_spin` uses:

| file | role |
| --- | --- |
| `prophet-v1-mag.pt` | energy model (`SpinCalculator`) |
| `prophet-v1-mag-head.pt` | magnetic-moment head (`MagmomPredictor`) |

They are **not** stored in git (137 MB each, see `.gitignore`).  Put them here
before running the prophet examples or building a binary that bundles them:

```bash
# from a local copy of the Kairos weights
cp /path/to/prophet-v1-mag.pt /path/to/prophet-v1-mag-head.pt models/

# or from the Hugging Face repository
python - <<'EOF'
from huggingface_hub import hf_hub_download
for name in ("prophet-v1-mag.pt", "prophet-v1-mag-head.pt"):
    hf_hub_download("kairosmaterials/prophet-v1-mag", name, local_dir="models")
EOF
```

`common.paths.model_file()` looks for them here (and inside the frozen bundle);
`PROPHET_MODEL` / `PROPHET_HEAD` override the paths, which is what a build made
with `--no-models` relies on.