# AmpliPiper Docker image

![Docker](https://github.com/nhmvienna/AmpliPiper/actions/workflows/docker-publish.yml/badge.svg)

## What is Docker?

[Docker](https://www.docker.com/) is a platform that help developers to get **faster**, more **generalized** and **packaged** deployments and distributions for their applications.

It can be regarded as a **DevOps** (Developement Operations) technology, as it is aimed at **easing and accelerating** shipping processes and times.

Docker is **cross-platform** (available on Linux, Windows and MacOs) and OS-independent as its main feature is allowing developers to build **isolated environments** (virtual machines) in which to run other softwares without having any conflict or versioning/platform-specific problem with the local machine.

It comes as a **command line** tool, but also, in its **Docker Desktop** implementation, as a Desktop application from which developers can comfortably check Docker's status on their computer.

We will talk about downloading and getting to run Docker in the next paragraph, but feel free to explore Docker resources [on their website](https://docs.docker.com/)!

## Get Docker

### Docker on Windows

[Docker on Windows](https://docs.docker.com/desktop/install/windows-install/) is only available as Docker Desktop, a Desktop application that manages _Docker engine_ (the technology that actually runs the virtual machines) through either **WSL2** or **Hyper-V**  as backend systems. 

You need Windows 10 or 11 to make Docker work and, if you're using a _Windows Subsystem Linux_ (WSL), make sure to have upgraded it to WSL2. 

For the installation, you should get the binary that best suits your Windows machine from the link at the beginning of this paragraph, and from there simply click on `Docker Desktop Installer.exe` and follow up as you are prompted by the installation interface.

### Docker on MacOS

[Docker on MacOS](https://docs.docker.com/desktop/install/mac-install/), as we said for Windows, can only be downloaded as Docker Desktop.

Docker, for now, supports only the three latest releases of macOS, and requires at least 4GB RAM to run.

You can get the binaries that most suit your platform from the link at the beginning of this paragraph, and then click on the `Docker.dmg` installer so obtained, dragging the Docker icon into the `Applications` folder.  Now you will have a `Docker.app`, which you will be able to initialize and set up following the prompts on the installation interface.

### Docker on Linux

There are [many supported platforms](https://docs.docker.com/engine/install/#supported-platforms) and many ways to install Docker for each one of them.

**1. Ubuntu and Debian-flavored Linux distributions**

You can simply run the following command:

```bash
sudo apt update && sudo apt install docker.io
```

And check the installation by running:

```bash
sudo docker --version
```

This method **will not** install `docker-compose` or `docker-buildx` and, for the sake of completeness, we will report here how to [get them](https://docs.docker.com/compose/install/linux/):

1. Start by getting Docker's official GPG Key:
```bash
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
```

2. Now add Docker and plugins to your `apt` repositories:
```bash
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
```

3. Now just install `compose` and `buildx`:
```bash
sudo apt-get install docker-compose-plugin docker-buildx-plugin
```

**2. CentOS-flavored Linux distributions**:

```bash
## Download rpm repo
sudo yum install -y yum-utils
sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo

## Install Docker engine
sudo yum install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

## Start Docker engine
sudo systemctl start docker
```

Other installations processes [here](https://docs.docker.com/engine/install/)

## Run AmpliPiper inside Docker

### 1. Use `docker compose` (RECOMMENDED)

Docker offers an integrated workflow that runs your images taking care of all the runtime arguments with a [`compose.yaml`](./compose.yaml) file. To run AmpliPiper container interactively, you just need to:

- Modify the  `USERDATA_PATH` variable in the [`.env`](./.env) to match with the portion of your file-system you want to inject into Docker

- > ⚠️: _If you are on Windows, you should **use PowerShell** and specify your paths in the `.env` file as `c:/path/to/your/data`._

- Run the following command **for Linux**:

    ```bash
    bash compose.sh
    ```

- Run the following command **for Windows PowerShell**:

    ```powershell
    .\compose.ps1
    ```

You'll find the volume you mounted in the container under `/app/userdata/`.


### 2. Use Docker image directly

For now AmpliPiper Docker image is available only for testing, and canbe obtained simply running from **PowerShell**/**CMD** (Windows), **Terminal** (macOS) and **integrated terminal** (Linux distros):

```bash
docker pull ghcr.io/nhmvienna/amplipiper:main
```

If it does not work on Linux for permission issues, you should run this as `sudo`.

You can now run the image interactively:

```bash
docker run -i \
    -t \
    ghcr.io/nhmvienna/amplipiper:main \
    /bin/bash
```

But the so-created container will not contain any data fromn your local file-system. 

In order to inject your local file-system into a Docker container you need to mount a volume (`-v` flag), using a mapping syntax that is very simple: `/your/local/path:/docker/container/path` (this is **only an example**! Replace `/your/local/path` and `/docker/container/path` with actual and valid paths).

To userdata AmpliPiper, you can clone this repository and then mount it inside the container as volume:

```bash
# Clone the repository
git clone https://github.com/nhmvienna/AmpliPiper.git

docker run -i \
    -t \
    -v ./AmpliPiper/:/app/userdata/ \
    ghcr.io/nhmvienna/amplipiper:main \
    /bin/bash
```

Now all the content in your local `AmpliPiper` folder is stored under `/app/userdata`.

## Test AmpliPiper

While you are inside the container, now, you can try the test already available in [AmpliPiper repo](https://github.com/nhmvienna/AmpliPiper) by running:

```bash
# define the path to your userdata folder
WD='/app/userdata'

# generate samples.csv file

## print header for samples.csv
printf "ID,PATH\n" >${WD}/testdata/data/samples.csv

## loop through input FASTQ files
for Filepath in ${WD}/testdata/reads/*fastq.gz; do

    ## get Filename
    Filename=${Filepath##*/}

    ## get File ID
    ID=${Filename%.fastq.gz*}

    ## print to samples.csv
    echo ${ID},${Filepath} >>${WD}/testdata/data/samples.csv
done

# run AmpliPiper
bash /app/shell/AmpliPiper.sh \
    -s /app/userdata/testdata/data/samples.csv \
    -p /app/userdata/testdata/data/primers.csv \
    -o /app/userdata/testdata/results/demo \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 200 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --force
```

You can find these instructions in [run.sh](./docker/run.sh).