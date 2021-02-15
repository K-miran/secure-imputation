echo "Cleaning data."
ls | grep -v clean.sh | grep -v setup.sh | xargs -Ifiles rm -f -r files 
