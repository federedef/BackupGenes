#!/usr/bin/env bash
. ~soft_bio_267/initializes/init_autoflow 

#Input variables.
exec_mode=$1 
add_opt=$2 

# Used Paths.
export input_path=`pwd`
export PATH=$input_path/scripts/aux_scripts:~soft_bio_267/programs/x86_64/scripts:$PATH
export autoflow_scripts=$input_path/scripts/autoflow_scripts
daemon_scripts=$input_path/scripts/daemon_scripts
export control_genes_folder=$input_path/control_genes
export output_folder=$SCRATCH/executions/BackupGenes
export output_folder_GraphPrioritizer=$SCRATCH/executions/GraphPrioritizer
report_folder=$output_folder/report

# Custom variables.
annotations="string_ppi disease phenotype biological_process string_ppi_textmining string_ppi_coexpression gene_hgncGroup"
kernels="ka rf ct el node2vec raw_sim"
integration_types="mean integration_mean_by_presence median max geometric_mean"
control_pos=$input_path'/control_pos'
control_neg=$input_path'/control_neg'

if [ "$exec_mode" == "download_translators" ] ; then
 
  
  ############################
  ## Obtain TRANSLATOR TABLES.
  mkdir -p ./translators

  # Downloading HGNC_symbol
  wget http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2022-04-01.txt -O ./translators/HGNC_symbol
  awk '{OFS="\t"}{print $2,$1}' ./translators/HGNC_symbol > ./translators/symbol_HGNC
  awk '{FS="\t";OFS="\t"}{print $19,$1}' ./translators/HGNC_symbol > ./translators/entrez_HGNC
  awk '{OFS="\t"}{print $2,$1}' ./translators/HGNC_symbol > ./translators/symbol_HGNC
 
elif [ "$exec_mode" == "control_preparation" ] ; then 

  $daemon_scripts/control_preparation.sh

elif [ "$exec_mode" == "control_type" ] ; then 

##################################################################
# OPTIONAL STAGE : SEE IF THE RELATION BACKUP-GENSEED IS SYMMETRIC
##################################################################
  . ~soft_bio_267/initializes/init_python
  filter_feature=$3 # Paralogs, Not_Paralogs, ".*"

  echo "$filter_feature"

  if [ $add_opt == "reverse" ] ; then 
      awk '{OFS="\t"}{print $2,$1,$3}' $control_genes_folder/backupgens/backup_gens | grep -w "$filter_feature" | cut -f 1,2 | aggregate_column_data -i - -x 1 -a 2 > ./control_pos
      awk '{OFS="\t"}{print $2,$1,$3}' $control_genes_folder/backupgens/non_backup_gens | grep -w "$filter_feature" | cut -f 1,2 | aggregate_column_data -i - -x 1 -a 2 > ./control_neg   
  elif [ $add_opt == "right" ] ; then 
      awk '{OFS="\t"}{print $1,$2,$3}' $control_genes_folder/backupgens/backup_gens | grep -w "$filter_feature" | cut -f 1,2 | aggregate_column_data -i - -x 1 -a 2 > ./control_pos
      awk '{OFS="\t"}{print $1,$2,$3}' $control_genes_folder/backupgens/non_backup_gens | grep -w "$filter_feature" | cut -f 1,2 | aggregate_column_data -i - -x 1 -a 2 > ./control_neg
  fi


elif [ "$exec_mode" == "ranking" ] ; then
  #########################################################
  # STAGE 2.2 OBTAIN RANKING FROM NON INTEGRATED KERNELS
  if [ -s $output_folder/rankings ] ; then
    rm -r $output_folder/rankings 
  fi
  mkdir -p $output_folder/rankings
  method=$2
  
  cat  $output_folder_GraphPrioritizer/similarity_kernels/*/*/ugot_path > $output_folder_GraphPrioritizer/similarity_kernels/ugot_path
  for annotation in $annotations ; do 
    for kernel in $kernels ; do 

      ugot_path="$output_folder_GraphPrioritizer/similarity_kernels/ugot_path"
      folder_kernel_path=`awk '{print $0,NR}' $ugot_path | sort -k 5 -r -u | grep "${annotation}_$kernel" | awk '{print $4}'`
      echo ${folder_kernel_path}
      if [ ! -z ${folder_kernel_path} ] ; then # This kernel for this annotation is done? 

        autoflow_vars=`echo " 
        \\$param1=$annotation,
        \\$kernel=$kernel,
        \\$input_path=$input_path,
        \\$folder_kernel_path=$folder_kernel_path,
        \\$input_name='kernel_matrix_bin',
        \\$production_seedgens=$production_seedgens,
        \\$control_pos=$control_pos,
        \\$control_neg=$control_neg,
        \\$output_name='non_integrated_rank',
        \\$method=$method,
        \\$geneseeds=$input_path/geneseeds
        " | tr -d [:space:]`

        AutoFlow -w $autoflow_scripts/ranking.af -V $autoflow_vars -o $output_folder/rankings/ranking_${kernel}_${annotation} -m 60gb -t 4-00:00:00 $3
      fi

    done
  done

elif [ "$exec_mode" == "integrated_ranking" ] ; then
  #########################################################
  # STAGE 2.4 OBTAIN RANKING FROM INTEGRATED KERNELS
  if [ -s $output_folder/integrated_rankings ] ; then
    rm -r $output_folder/integrated_rankings # To not mix executions.
  fi
  mkdir -p $output_folder/integrated_rankings
  cat  $output_folder_GraphPrioritizer/integrations/*/*/ugot_path > $output_folder_GraphPrioritizer/integrations/ugot_path # What I got?

  method=$2

  for integration_type in ${integration_types} ; do 
    for kernel in $kernels ; do 

      ugot_path="$output_folder_GraphPrioritizer/integrations/ugot_path"
      folder_kernel_path=`awk '{print $0,NR}' $ugot_path | sort -k 5 -r -u | grep -w "${integration_type}_$kernel" | awk '{print $4}'`
      echo ${folder_kernel_path}
      if [ ! -z ${folder_kernel_path} ] ; then # This kernel for this integration_type is done? 

        autoflow_vars=`echo " 
        \\$param1=$integration_type,
        \\$kernel=$kernel,
        \\$input_path=$input_path,
        \\$folder_kernel_path=$folder_kernel_path,
        \\$input_name='general_matrix',
        \\$production_seedgens=$production_seedgens,
        \\$control_pos=$control_pos,
        \\$control_neg=$control_neg,
        \\$output_name='integrated_rank',
        \\$method=$method,
        \\$geneseeds=$input_path/geneseeds
        " | tr -d [:space:]`

        AutoFlow -w $autoflow_scripts/ranking.af -V $autoflow_vars -o $output_folder/integrated_rankings/ranking_${kernel}_${integration_type} -m 60gb -t 4-00:00:00 $3
      fi

    done
  done

#########################################################
# STAGE 3 OBTAIN REPORT FROM RESULTS
#########################################################

elif [ "$exec_mode" == "report" ] ; then 
  source ~soft_bio_267/initializes/init_python
  source ~soft_bio_267/initializes/init_R
  html_name=$2
  check=$3
  
  # #################################
  # Setting up the report section #
  find $report_folder/ -mindepth 2 -delete
  find $output_folder/ -maxdepth 1 -type f -delete

  mkdir -p $report_folder/ranking_report
  mkdir -p $report_folder/img

  declare -A original_folders

  original_folders[non_integrated_rank_summary]='rankings'
  original_folders[non_integrated_rank_measures]='rankings'
  original_folders[non_integrated_rank_cdf]='rankings'
  original_folders[non_integrated_rank_pos_cov]='rankings'
  original_folders[non_integrated_rank_positive_stats]='rankings'

  original_folders[integrated_rank_summary]='integrated_rankings'
  original_folders[integrated_rank_measures]='integrated_rankings'
  original_folders[integrated_rank_cdf]='integrated_rankings'
  original_folders[integrated_rank_pos_cov]='integrated_rankings'
  original_folders[integrated_rank_positive_stats]='integrated_rankings'
  
  # Here the data is collected from executed folders.
  for file in "${!original_folders[@]}" ; do
    original_folder=${original_folders[$file]}
    count=`find $output_folder/$original_folder -maxdepth 3 -mindepth 3 -name $file | wc -l`
    if [ "$count" -gt "0" ] ; then
      echo "$file"
      cat $output_folder/$original_folder/*/*/$file > $output_folder/$file
    fi
  done 

  ##########################
  # Processing all metrics #
  declare -A references

  references[non_integrated_rank_summary]='Sample,Net,Kernel'
  references[non_integrated_rank_pos_cov]='Sample,Net,Kernel'
  references[non_integrated_rank_positive_stats]='Sample,Net,Kernel,group_seed'

  references[integrated_rank_summary]='Sample,Integration,Kernel'
  references[integrated_rank_pos_cov]='Sample,Integration,Kernel'
  references[integrated_rank_positive_stats]='Sample,Integration,Kernel,group_seed'

  references[annotation_grade_metrics]='Gene_seed'

  for metric in non_integrated_rank_summary integrated_rank_summary non_integrated_rank_pos_cov integrated_rank_pos_cov non_integrated_rank_positive_stats integrated_rank_positive_stats ; do
    if [ -s $output_folder/$metric ] ; then
      echo "$output_folder/$metric"
      create_metric_table $output_folder/$metric ${references[$metric]} $report_folder/ranking_report/parsed_${metric} 
    fi
  done

  if [ -s $output_folder/non_integrated_rank_measures ] ; then
     echo -e "annot_kernel\tannot\tkernel\trank\tacc\ttpr\tfpr\tprec\trec" | \
     cat - $output_folder/non_integrated_rank_measures > $report_folder/ranking_report/non_integrated_rank_measures
  fi

    if [ -s $output_folder/integrated_rank_measures ] ; then
    echo -e "integration_kernel\tintegration\tkernel\trank\tacc\ttpr\tfpr\tprec\trec" | \
     cat - $output_folder/integrated_rank_measures > $report_folder/ranking_report/integrated_rank_measures
  fi

  if [ -s $output_folder/non_integrated_rank_cdf ] ; then
     echo -e "annot_kernel\tannot\tkernel\tcandidate\tscore\trank\tcummulative_frec\tgroup_seed"| \
     cat - $output_folder/non_integrated_rank_cdf > $report_folder/ranking_report/non_integrated_rank_cdf
  fi

  if [ -s $output_folder/integrated_rank_cdf ] ; then
     echo -e "integration_kernel\tintegration\tkernel\tcandidate\tscore\trank\tcummulative_frec\tgroup_seed"| \
     cat - $output_folder/integrated_rank_cdf > $report_folder/ranking_report/integrated_rank_cdf
  fi
  
  if [ -z "$check" ] ; then
    echo "---------------------------------------"
    echo " Now it is necessary some information of the process "
    echo "data version?"
    read version
    echo "extra info?"
    read extra_info
    name_dir=`date +%d_%m_%Y`
    mkdir ./report/HTMLs/$name_dir
    # Create preprocess file
    echo -e " Data version:\t$version\nExtra Info:\t$extra_info " > ./report/HTMLs/$name_dir/info_preprocess
  else 
    echo "Reports to check available"
  fi
 ###################
  # Obtaining HTMLS #
  report_html -t ./report/templates/ranking_report.mko -d `ls $report_folder/ranking_report/* | tr -s [:space:] "," | sed 's/,*$//g'` -o "report_algQuality$html_name"

  if [ -z "$check" ] ; then
    mv ./report_algQuality$html_name.html ./report/HTMLs/$name_dir/
  fi

#########################################################
# STAGE TO CHECK AUTOFLOW IS RIGHT
#########################################################
elif [ "$exec_mode" == "check" ] ; then
  #STAGE 3 CHECK EXECUTION
  for folder in `ls $output_folder/$add_opt/` ; do 
    if [ -d $output_folder/$add_opt/$folder ] ; then
      echo "$folder"
      flow_logger -w -e $output_folder/$add_opt/$folder -r all
    fi
  done  
fi
