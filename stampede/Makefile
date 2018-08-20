APP=mattmiller899-demultiplex-app-0.0.1

clean:
	find . \( -name \*.out -o -name \*.param \) -exec rm {} \;

up:
	copy_from_manifest.py

template:
	jobs-template $(APP) > job.json

job:
	jobs-submit -F job.json

app:
	apps-addupdate -F app.json

test:
	sbatch test.sh
