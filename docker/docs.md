# Technical documentation for AmpliPiper Docker image

The right context to build AmpliPiper Docker image can be found [here](https://github.com/AstraBert/amplipiper)

## Dockerfile

Based on [`condaforge/miniforge3`](https://hub.docker.com/r/mambaorg/micromamba) Docker image, the Dockerfile includes a:

- Single-stage build from the aforementioned image
- Working directory is set as `/app/`
- Context transferring from local filesystem to Docker image
- Execution of the setup script (find it at [`shell/setup.sh`](./shell/setup.sh))

## CI/CD

Continuous-Integration/Continous-Developement is maintained through [the base repo for this Docker image](https://github.com/AstraBert/amplipiper) and managed with [**GitHub Container Registry**](https://ghcr.io/). the image will be availble under: `ghcr.io/astrabert/amplipiper` and we advise you pull the `main` tag, as it is always on track with the modifications operated on the `main` branch of the repository (the one on which the image is based).

## Build locally

Execute the code to build the image locally:

```bash
git clone https://github.com/AstraBert/amplipiper.git
cd amplipiper
docker build . -t YOUR-USER-NAME/IMAGE-NAME:tag
```