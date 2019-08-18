#!/usr/bin/env bash

target_dir=/lustre/scratch117/cellgen/cellgeni/forecast;
file_map_name="file_map.txt"
mkdir -p ${target_dir};
counter=0

copy_file(){
    counter=$((${counter}+1));
    cp -n $1 ${target_dir}/input.${counter};
    # echo $1;
}

create_file_mapping(){
    counter=$((${counter}+1));
    echo $1 ${target_dir}/input.${counter} | tee -a ${file_map_name};
}

run_wge_to_ccid_on_file(){
    source activate forecast
    bsub -q normal -n 1 -R "span[hosts=1]" -M 20 -R "select[mem>=20] rusage[mem=20]" ./wge_to_ccid.py
}

parse_ccds(){
    for file in $1/*
    do
        if [[ ${file} == *reads.txt ]]
        then
#            copy_file ${file};
            create_file_mapping ${file};
        fi
    done
}

parse_species_dir(){
    for dir in $1/*;
    do
        if [[ ${dir} == *CCDS* ]]
        then
            parse_ccds ${dir};
        fi
    done
}

main(){
    root_directory="/lustre/scratch117/cellgen/team227/FORECasT_profiles_for_AK"
    human_dir="${root_directory}/human"
    mouse_dir="${root_directory}/mouse"

    parse_species_dir ${human_dir}
    parse_species_dir ${mouse_dir}
    echo ${counter}
}

main
