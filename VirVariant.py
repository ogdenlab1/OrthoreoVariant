import pandas as pd
import os
import fnmatch
import argparse
from pathlib import Path

#to do: add ability to parse viral annotation input for gene-level annotation

def getCoverage(cov_file, sample, df):
    df_depth = pd.read_csv(cov_file, sep="\t", header=0)
    total_depth = sum(df_depth["Coverage"])
    df.loc[df['sample'] == sample, ['total_nts']] = total_depth
    return df

def getVCF(vcf_file):
    vcf = pd.read_csv(vcf_file, skiprows=18, sep="\t", header=0, index_col=False, names=["genome",
                                                                                          "position",
                                                                                          "ID",
                                                                                          "reference",
                                                                                          "variant",
                                                                                          "qual",
                                                                                          "filter",
                                                                                          "info"])
    vcf[['raw_depth',
         'frequency',
         'strand_bias',
         'DP4']] = vcf['info'].apply(lambda x: pd.Series(x.split(';')))
    vcf['raw_depth'] = vcf['raw_depth'].str[3:]
    vcf['frequency'] = vcf['frequency'].str[3:]
    vcf['DP4'] = vcf['DP4'].str[4:]
    vcf[['ref_f_count', 'ref_r_count', 'variant_f_count', 'variant_r_count']] = vcf['DP4'].apply(
        lambda x: pd.Series(x.split(',')))
    vcf = vcf.drop(columns=['ID', "qual", "filter", 'info', 'strand_bias', 'DP4'])
    vcf = vcf[['genome', 'position', 'reference', 'variant', 'frequency', 'raw_depth', 'ref_f_count', 'ref_r_count',
               'variant_f_count', 'variant_r_count']]
    vcf['position'] = pd.to_numeric(vcf['position'])
    vcf['frequency'] = pd.to_numeric(vcf['frequency'])
    vcf['raw_depth'] = pd.to_numeric(vcf['raw_depth'])
    vcf['ref_f_count'] = pd.to_numeric(vcf['ref_f_count'])
    vcf['ref_r_count'] = pd.to_numeric(vcf['ref_r_count'])
    vcf['variant_f_count'] = pd.to_numeric(vcf['variant_f_count'])
    vcf['variant_r_count'] = pd.to_numeric(vcf['variant_r_count'])
    vcf['variant_total'] = vcf['variant_f_count'] + vcf['variant_r_count']
    return vcf


def get_variant_type(reference, variant):
    type = ""
    if reference == "A" and variant == "G":
        type = "transition"
    elif reference == "G" and variant == "A":
        type = "transition"
    elif reference == "C" and variant == "T":
        type = "transition"
    elif reference == "T" and variant == "C":
        type = "transition"
    else:
        type = "transversion"
    return type


def calculateVariants(df, out_report, sample):
    variant_nts = df['variant_total'].sum()
    transition_nts = df.loc[df['variant_type'] == "transition", 'variant_total'].sum()
    transversion_nts = df.loc[df['variant_type'] == "transversion", 'variant_total'].sum()
    AtoG_nts = df.loc[((df['reference'] == "A") & (df['variant'] == "G")), 'variant_total'].sum()
    GtoA_nts = df.loc[((df['reference'] == "G") & (df['variant'] == "A")), 'variant_total'].sum()
    CtoT_nts = df.loc[((df['reference'] == "C") & (df['variant'] == "T")), 'variant_total'].sum()
    TtoC_nts = df.loc[((df['reference'] == "T") & (df['variant'] == "C")), 'variant_total'].sum()
    AtoC_nts = df.loc[((df['reference'] == "A") & (df['variant'] == "C")), 'variant_total'].sum()
    CtoA_nts = df.loc[((df['reference'] == "C") & (df['variant'] == "A")), 'variant_total'].sum()
    AtoT_nts = df.loc[((df['reference'] == "A") & (df['variant'] == "T")), 'variant_total'].sum()
    TtoA_nts = df.loc[((df['reference'] == "T") & (df['variant'] == "A")), 'variant_total'].sum()
    CtoG_nts = df.loc[((df['reference'] == "C") & (df['variant'] == "G")), 'variant_total'].sum()
    GtoC_nts = df.loc[((df['reference'] == "G") & (df['variant'] == "C")), 'variant_total'].sum()
    GtoT_nts = df.loc[((df['reference'] == "G") & (df['variant'] == "T")), 'variant_total'].sum()
    TtoG_nts = df.loc[((df['reference'] == "T") & (df['variant'] == "G")), 'variant_total'].sum()
    unique_variants = len(df)
    out_report.loc[out_report['sample'] == sample, ['unique_variants']] = unique_variants
    out_report.loc[out_report['sample'] == sample, ['variant_nts']] = variant_nts
    out_report.loc[out_report['sample'] == sample, ['transition_nts']] = transition_nts
    out_report.loc[out_report['sample'] == sample, ['transversion_nts']] = transversion_nts
    out_report.loc[out_report['sample'] == sample, ['AtoG_nts']] = AtoG_nts
    out_report.loc[out_report['sample'] == sample, ['GtoA_nts']] = GtoA_nts
    out_report.loc[out_report['sample'] == sample, ['CtoT_nts']] = CtoT_nts
    out_report.loc[out_report['sample'] == sample, ['TtoC_nts']] = TtoC_nts
    out_report.loc[out_report['sample'] == sample, ['AtoC_nts']] = AtoC_nts
    out_report.loc[out_report['sample'] == sample, ['CtoA_nts']] = CtoA_nts
    out_report.loc[out_report['sample'] == sample, ['AtoT_nts']] = AtoT_nts
    out_report.loc[out_report['sample'] == sample, ['TtoA_nts']] = TtoA_nts
    out_report.loc[out_report['sample'] == sample, ['CtoG_nts']] = CtoG_nts
    out_report.loc[out_report['sample'] == sample, ['GtoC_nts']] = GtoC_nts
    out_report.loc[out_report['sample'] == sample, ['GtoT_nts']] = GtoT_nts
    out_report.loc[out_report['sample'] == sample, ['TtoG_nts']] = TtoG_nts
    return out_report, df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("Sample_List",
                        help="A tab delineated file with each line containing sample names. Last line is empty.")
    parser.add_argument("Working_Directory", help="Path to directory containing data to align.")
    parser.add_argument("Experiment", help="Experiment name.")
    parser.add_argument("--freq", help="Variant frequency cutoff for filtering. Decimal between 0 and 1.", default=0, required=False)
    parser.add_argument("--file_tag", help="File naming tag can denote filters used for variant isolation.", default="", required=False)
    args = parser.parse_args()
    sample_list = [line.rstrip('\n') for line in open(str(args.Sample_List))]
    wd = Path(args.Working_Directory)
    exp = str(args.Experiment)
    if args.file_tag:
        tag = "_" + str(args.file_tag)
    else:
        tag = ""
    if args.freq:
        freq_cutoff = float(args.freq)
    else:
        freq_cutoff = 0
    report = pd.DataFrame(columns=['sample',
                                   'unique_variants',
                                   'variant_nts',
                                   "total_nts",
                                   "transition_nts",
                                   "transversion_nts",
                                   "AtoG_nts",
                                   "GtoA_nts",
                                   "CtoT_nts",
                                   "TtoC_nts",
                                   "AtoT_nts",
                                   "TtoA_nts",
                                   "AtoC_nts",
                                   "CtoA_nts",
                                   "CtoG_nts",
                                   "GtoC_nts",
                                   "GtoT_nts",
                                   "TtoG_nts",
                                   "mutation_freq",
                                   "transition_freq",
                                   "transversion_freq",
                                   "AtoG_freq",
                                   "GtoA_freq",
                                   "CtoT_freq",
                                   "TtoC_freq",
                                   "AtoT_freq",
                                   "TtoA_freq",
                                   "AtoC_freq",
                                   "CtoA_freq",
                                   "CtoG_freq",
                                   "GtoC_freq",
                                   "GtoT_freq",
                                   "TtoG_freq"
                                   ])
    report['sample'] = sample_list
    for sample in sample_list:
        report = getCoverage(Path(wd, f"{sample}_bowtie2_coverage.txt"), sample, report)
        sample_vcf = getVCF(Path(wd, f"{sample}.vcf"))
        sample_vcf_filter = sample_vcf[sample_vcf['frequency'] >= freq_cutoff]
        sample_vcf_filter['variant_type'] = sample_vcf_filter[['reference', 'variant']].apply(lambda x: get_variant_type(*x), axis=1)
        report, variant_output = calculateVariants(sample_vcf_filter, report, sample)
        report['mutation_freq'] = (report['variant_nts'] / report['total_nts'])
        report['transition_freq'] = report['transition_nts'] / report['total_nts']
        report['transversion_freq'] = report['transversion_nts'] / report['total_nts']
        report['AtoG_freq'] = report['AtoG_nts'] / report['total_nts']
        report['GtoA_freq'] = report['GtoA_nts'] / report['total_nts']
        report['AtoC_freq'] = report['AtoC_nts'] / report['total_nts']
        report['CtoA_freq'] = report['CtoA_nts'] / report['total_nts']
        report['AtoT_freq'] = report['AtoT_nts'] / report['total_nts']
        report['TtoA_freq'] = report['TtoA_nts'] / report['total_nts']
        report['CtoG_freq'] = report['CtoG_nts'] / report['total_nts']
        report['GtoC_freq'] = report['GtoC_nts'] / report['total_nts']
        report['CtoT_freq'] = report['CtoT_nts'] / report['total_nts']
        report['TtoC_freq'] = report['TtoC_nts'] / report['total_nts']
        report['GtoT_freq'] = report['GtoT_nts'] / report['total_nts']
        report['TtoG_freq'] = report['TtoG_nts'] / report['total_nts']
        #save outputs
        variant_output.to_csv(Path(wd, f"{sample}{tag}_variants.txt"), sep="\t", index=False)
    report.to_csv(Path(wd, f"{exp}_variant_summary.txt"), sep="\t", index=False)


if __name__ == '__main__':
    main()
