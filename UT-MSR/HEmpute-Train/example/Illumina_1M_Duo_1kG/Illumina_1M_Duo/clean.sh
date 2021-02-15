echo "Cleaning data."
ls | grep -v clean.sh | grep -v setup_run_models.sh | xargs -Ifiles rm -f -r files 
