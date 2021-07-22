source ~/.bashrc
      cd ~/projs/scATACseq_MLL

    for num in num102 num104 num105 num108 num112 num113 num114 num115 numC1 numC2 numC3 numC4 numC5
    do
    snaptools snap-del \
        --snap-file snaps/${num}.snap \
        --session-name PM
    done
echo "snap-del is done."

      for num in num102 num104 num105 num108 num112 num113 num114 num115 numC1 numC2 numC3 numC4 numC5
      do
      snaptools snap-add-pmat \
    	  --snap-file snaps/${num}.snap \
    	  --peak-file peaks.combined.bed
    	done
echo "snap-add-pmat is done."
