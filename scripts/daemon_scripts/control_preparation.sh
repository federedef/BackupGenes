source ~soft_bio_267/initializes/init_R
source ~soft_bio_267/initializes/init_python

  # Backup Controls #
  ###################
  mkdir -p $control_genes_folder/backupgens/processed_data 

  # Process data 
  
  # POSITIVE CONTROL #
  echo "POSITIVE CONTROL"
  # AdHoc added backups
  cp $control_genes_folder/backupgens/data/AdHoc_Backups $control_genes_folder/backupgens/processed_data/AdHoc_Backups
  # Big Papi.
  grep -w 'All 6\|no *' $control_genes_folder/backupgens/data/Big_Papi | awk '{FS="\t";OFS="\t"}{if ( $6 < 0.05 ) print $1,$2}' | sort -k1 | uniq | grep -v -e '^$' > $control_genes_folder/backupgens/data/filtered_Big_Papi
  # Digenic Paralog.
  awk '{FS="\t"}{if ( $2 <= 0.05 && $3 <= 0.05 && $4 <= 0.05 && $5 <= 0.05 && $6 <= 0.05 && $7 <= 0.05 && $8 <= 0.05 && $9 <= 0.05 && $10 <= 0.05 && $11 <= 0.05 && $12 <= 0.05) print $1}' $control_genes_folder/backupgens/data/Digenic_Paralog | \
   tr -s ";" "\t" | tr -d "\"" | grep -v -e '^$' > $control_genes_folder/backupgens/data/filtered_Digenic_Paralog

  standard_name_replacer -u -I ./translators/symbol_HGNC -i $control_genes_folder/backupgens/data/filtered_Big_Papi -c 1,2 > $control_genes_folder/backupgens/processed_data/Big_Papi
  standard_name_replacer -u -I ./translators/symbol_HGNC -i $control_genes_folder/backupgens/data/filtered_Digenic_Paralog -c 1,2 > $control_genes_folder/backupgens/processed_data/Digenic_Paralog
  
  cat $control_genes_folder/backupgens/processed_data/* | sort | uniq -u  > $control_genes_folder/backupgens/backup_gens

  # NEGATIVE CONTROL #
  echo "NEGATIVE CONTROL"
  #grep -w 'All 6' $control_genes_folder/backupgens/data/Big_Papi | awk '{FS="\t";OFS="\t"}{if ( $7 <= 0.05) print $1,$2}' > $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control
  grep -w 'All 6\|no *' $control_genes_folder/backupgens/data/Big_Papi|  awk '{FS="\t";OFS="\t"}{if ( $6 >= 0.5 ) print $1,$2}' | sort -k1 | uniq -c | grep -w "6" | sed "s/6//1" | tr -d " "  > $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control
  awk '{FS="\t"}{if ( $2 >= 0.5 && $3 >= 0.5 && $4 >= 0.5 && $5 >= 0.5 && $6 >= 0.5 && $7 >= 0.5 && $8 >= 0.5 && $9 >= 0.5 && $10 >= 0.5 && $11 >= 0.5 && $12 >= 0.5) print $1}' $control_genes_folder/backupgens/data/Digenic_Paralog | \
   tr -s ";" "\t" | tr -d "\"" >> $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control
  grep -v -e '^$' $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control > tmp
  mv tmp $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control
  standard_name_replacer -u -I ./translators/symbol_HGNC -i $control_genes_folder/backupgens/data/filtered_Big_Papi_negative_control -c 1,2 | awk '{if (!( $1 == $2 )) print $0 }' > $control_genes_folder/backupgens/non_backup_gens
  

  echo "Obtaining Paralogs genes"
  # Finally add new column indicating which pairs are paralogs in NEGATIVE AND POSITIVE CONTROLS.
  which_are_paralogs.R -i $control_genes_folder/backupgens/backup_gens -o "$control_genes_folder/backupgens" -O "backup_gens"
  echo "positive tagged"
  which_are_paralogs.R -i $control_genes_folder/backupgens/non_backup_gens -o "$control_genes_folder/backupgens" -O "non_backup_gens"

