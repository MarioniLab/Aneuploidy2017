#!/usr/bin/bash
if [ ! -e "PaperSupplementary.html" ]
	then
		echo "There is no PaperSupplementary.html file here - are you running the script from the right place?"
		echo "If you still want to download here, run 'touch aneu_report.html' and then this script again."
	else


		echo "Getting processed data"
		curl -k https://jmlab-gitlab.cri.camres.org/Jonny/Aneuploidy2017Data/raw/master/aneu_proc.tar.gz > aneu_proc.tar.gz
		echo "Unzipping processed data"
		tar -xzf aneu_proc.tar.gz

		echo "Getting raw data"
		curl -k https://jmlab-gitlab.cri.camres.org/Jonny/Aneuploidy2017Data/raw/master/aneu_raw.tar.gz > aneu_raw.tar.gz
		echo "Unzipping raw data"
		tar -xzf aneu_raw.tar.gz

fi
