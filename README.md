
## Instructions 

Please see guidelines below to run this code.

### Setup conda environment 

1. For more information on how to create and activate conda activate see this link https://www.anaconda.com/docs/getting-started/working-with-conda/environments. However, these are the two basic commands for creating and activating conda environments.  

        conda create -n <ENV_NAME>

        conda activate <ENV_NAME>

    NB: Replace  <ENV_NAME>  with the name of your environment


2. Use `which python` and `which pip` to make sure the environment is using the appropriate python


3. Next, download the git repo - ideally in a new folder. 


4. Make sure you're in the same directory where the `setup.py` is located before running the following command:

        pip3 install -e . 

    The command above downloads relevant dependencies to run `probeSeq` module

5. To see which libraries/packages have been download use the following command:
         
        pip3 list

### How to run / open juptyer notebook

6. Change your directory to where `processGenbank.ipynb` is located on your computer, and then use 
`juptyter notebook` into your terminal. This will trigger your browser and open the probePanel directory. 


7. Click on `processGenbank.ipynb` to open the notebook.