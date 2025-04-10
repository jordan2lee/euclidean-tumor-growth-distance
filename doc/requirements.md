# Requirements

## 1. Required Applications

+ Install [Python 3+](https://www.python.org/downloads/)
+ Install [Docker Desktop](https://docs.docker.com/get-started/get-docker/) (or Docker)
+ Install [GDC Data Transfer Tool Client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)

    Get app to download GDC data, download the GDC Data Transfer Tool Client for your system at the above link. Below is an example for Ubuntu.
    ```bash
    curl -O https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
    unzip gdc-client_v1.6.1_Ubuntu_x64.zip
    ```

## 2. Create an environment and install dependencies
```bash
python3 -m venv venv; source venv/bin/activate
```

```bash
pip install -r requirements.txt 
```

## 3. Sign into Applications
Docker Sign In (if not already)

```bash
docker login
```

## 4. GDC Client Token
Token setup. Go to https://portal.awg.gdc.cancer.gov/ get token under profile and save as `secrets/gdc-awg.token`

## 4. Download Required Data
