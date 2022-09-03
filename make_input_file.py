"""
Creates an input file from the directory (python 2.7)
"""

from os import listdir

def getSampleObj():
    return {
        "indiv_id": "",
        "sample_id": "",
        "lib_id": "",
        "rg_id": "",
        "plat_unit": "ILLUMINA",
        "platform": "",
        "plat_model": "",
        "run_date": "05/18/2022",
        "center": "AZENTA",
        "r1": "",
        "r2": ""
    }

if __name__ == "__main__":
    DATA_DIR = "/restricted/projectnb/pcga/JnJ_invitro/Data/" # actual data
    # DATA_DIR = "/Users/danyfu/Dropbox/CompBioMed/nextflow/RNA_Seq/JnJInvitro" # for testing
    TAB = "\t"
    FILES = [f for f in listdir(DATA_DIR) if f.endswith(".gz")]
    all_samples = {}
    for file_name in FILES: # HBEC-Control-1_R1_001.fastq.gz
        sample_name = file_name.split("_")[0] # HBEC-Control-1
        indiv_name = sample_name[:-2] # HBEC-Control
        if sample_name not in all_samples:
            all_samples[sample_name] = getSampleObj()
            all_samples[sample_name]["indiv_id"] = indiv_name
            all_samples[sample_name]["sample_id"] = sample_name
        if "_R1_" in file_name:
            all_samples[sample_name]["r1"] = "/restricted/projectnb/pcga/JnJ_invitro/Data/"+file_name
        elif "_R2_" in file_name:
            all_samples[sample_name]["r2"] = "/restricted/projectnb/pcga/JnJ_invitro/Data/"+file_name

    input_file = open('jnj_invitro_input.txt', 'w')
    input_file.writelines('INDIVIDUAL_ID{tab}SAMPLE_ID{tab}LIBRARY_ID{tab}RG_ID{tab}PLATFORM_UNIT{tab}PLATFORM{tab}PLATFORM_MODEL{tab}RUN_DATE{tab}CENTER{tab}R1{tab}R2'.format(tab=TAB))
    input_file.writelines('\n')

    for k, v in all_samples.iteritems():
        input_file.writelines('{indiv_id}{tab}{sample_id}{tab}{lib_id}{tab}{rg_id}{tab}{plat_unit}{tab}{platform}{tab}{plat_model}{tab}{run_date}{tab}{center}{tab}{r1}{tab}{r2}'
                              .format(tab=TAB, indiv_id=v["indiv_id"], sample_id=v["sample_id"], lib_id=v["lib_id"], rg_id=v["rg_id"],
                                      plat_unit=v["plat_unit"], platform=v["platform"], plat_model=v["plat_model"],
                                      run_date=v["run_date"], center=v["center"], r1=v["r1"], r2=v["r2"]))
        input_file.writelines('\n')
