#!/usr/bin/bash
if [ ! -e "aneu_report.html" ]
	then
		echo "There is no aneu_report.html file here - are you running the script from the right place?"
		echo "If you still want to download here, run 'touch aneu_report.html' and then this script again."
	else
		if [ ! -d "raw_data" ]
			then 
			mkdir raw_data
		fi

		if [ ! -d "proc_data" ]
			then 
			mkdir proc_data
		fi

		cd raw_data
		echo "Beginning ftp"
		wget --user=jmlabftp --password='HOBICAmeer6' ftp://ftp2.cruk.cam.ac.uk/aneuploidy/aneu_raw.tar.gz
		echo "Unzipping files"
		tar -xzf aneu_raw.tar.gz

fi
