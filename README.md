# PredAPP source code
Web server: http://PredAPP.xialab.info

Email: [jfxia@ahu.edu.cn](mailto:jfxia@ahu.edu.cn)



The repo is organised as follows:

- `/static`: Folder containing all static files in the project

- `/Models`: folder containing the 9 machine-learning models employed for the predictions in `.pkl` format.

- `/templates`: folder with the HTML files of the project
- `/datasets`: 
  - `Training.fasta`: Training dataset for PredAPP in Fasta format
  - `Test.fasta`: Independent test dataset for PredAPP in Fasta format
- `requirements.txt`: environment file with all dependencies needed to run the project

- `app.py`: Flask application for PredAPP

- `Features.py`: Python script for feature engineering
- `PredAPP.py`: Python script for machine learning model construction and evaluation



## Installation
- Requirement
  
  OS：
  
  - `Windows` ：Windows7 or later
  
  - `Linux`：Ubuntu 16.04 LTS or later
  
  Python：
  
  - `Python` >= 3.6

- Download `PredAPP`to your computer

  ```bash
  git clone https://github.com/xialab-ahu/PredAPP.git
  ```

- open the dir and install `requirement.txt` with `pip`

  ```
  cd PredAPP
  pip install -r requirement.txt
  ```

## Run

```shell
python -m flask run
```

Open the browser, open the url：http://127.0.0.1/5000

