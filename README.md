# Top-down Organization

This repository contains scripts and Docker setups for neuroscience data analysis.

## Repository layout

- **Containers/** &ndash; Dockerfiles and `docker-compose` configurations for
  running Jupyter Lab with OpenCV and an RStudio environment. Customize the port
  numbers (`XXXX`) and mount points to suit your machine before launching the
  containers.
- **Patch_clamp_analysis/** &ndash; Jupyter notebooks and R scripts for analyzing
  patch clamp experiments. They reference data under `./Data/` so place the
  necessary recordings there or update the paths inside the notebooks/scripts.
- **Projection_analysis/** &ndash; Materials for analyzing histological data,
  including notebooks for intensity mapping and R scripts for statistics.
- **README.md** &ndash; You are reading it. Brief usage information and an
  overview of the repository.
- **LICENSE** &ndash; MIT License.

## Getting started

1. Clone this repository and prepare your data directory (e.g. `./Data`).
2. Adjust the `docker-compose.yml` files in `Containers/OpenCV` and
   `Containers/Rstudio` &ndash; set port numbers and bind mounts.
3. Launch the environment:

   ```bash
   cd Containers/OpenCV
   docker compose up -d
   ```

   Jupyter Lab will be available on the port you specified. The RStudio
   environment can be started in the same way from `Containers/Rstudio`.
4. Open the notebooks or R scripts under `Patch_clamp_analysis` or
   `Projection_analysis` and run the analyses.

## Notes

- The repository does not contain raw data. Adjust paths to your local dataset
  as needed.
- There is currently no automated testing or CI setup; scripts are run
  interactively in the Docker containers.

Feel free to extend this README with additional instructions as new analyses are
added.
## Tips for new contributors
- Review the structure under `Patch_clamp_analysis` to see how patch clamp logs and ABF files are loaded.
- Explore `Projection_analysis` notebooks for examples of image-based analyses and statistics.
- Use the Docker containers to ensure dependencies (Python, R, CUDA) match across machines.
- No automated tests are provided, so execute notebooks or scripts directly to confirm they run in your setup.

