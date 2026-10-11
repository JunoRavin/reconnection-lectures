# Magnetic Reconnection with Gkeyll

Two computational lectures on magnetic reconnection: Sweet–Parker scaling and
fast reconnection, worked through with multi-fluid simulations in
[Gkeyll](https://github.com/gkeyllorg/gkeyll) and analyzed with
[postgkyl](https://github.com/gkeyllorg/postgkyl).

Read the lectures at <https://junoravin.github.io/reconnection-lectures>.

## Layout

- `notebooks/`: the lecture notebooks ([marimo](https://marimo.io)), the source of the lectures
- `lectures/`: the same notebooks exported with outputs (`.ipynb`); the website is built from these
- `inputs/`: Gkeyll input files for the simulations the lectures analyze
- `index.qmd`, `install.qmd`, `running.qmd`: the remaining chapters
- `scripts/export-lectures.sh`: re-runs the notebooks and regenerates the `.ipynb` files

## Updating the lectures

Edit a notebook with `marimo edit notebooks/<name>.py`, then regenerate its outputs
(this needs the simulation data):

```bash
GEM_DATA_DIR=/path/to/runs scripts/export-lectures.sh
quarto preview
```

Pushing to `main` rebuilds and deploys the site from the committed `.ipynb` files.
