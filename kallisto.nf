/*
 * Copyright (c) 2015-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Kallisto-NF'.
 *
 *   Kallisto-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Kallisto-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Kallisto-NF.  If not, see <http://www.gnu.org/licenses/>.

     Usage:
     nextflow run kallisto.nf \
     --reads s3://lifebit-featured-datasets/pipelines/kallisto-nf-data/reads \
     --transcriptome s3://lifebit-featured-datasets/pipelines/kallisto-nf-data/transcriptome/transcriptome.fa \
     --experiment s3://lifebit-featured-datasets/pipelines/kallisto-nf-data/experiment/hiseq_info.txt \
     -with-docker -resume
 */

/*
 * Main Kallisto-NF pipeline script
 *
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com>
 */


params.transcriptome = "$baseDir/tutorial/transcriptome/transcriptome.fa"
params.name          = "RNA-Seq Abundance Analysis"
params.reads         = "$baseDir/tutorial/reads"
params.readsExtension="fastq"
allReads="${params.reads}/*.${params.readsExtension}"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.experiment    = "$baseDir/tutorial/experiment/hiseq_info.txt"
params.output        = "results/"
params.multiqc_config = "$baseDir/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)
params.skip_multiqc = false


log.info "K A L L I S T O - N F  ~  version 0.9"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "reads                  : ${params.reads}"
log.info "transcriptome          : ${params.transcriptome}"
log.info "fragment length        : ${params.fragment_len} nt"
log.info "fragment SD            : ${params.fragment_sd} nt"
log.info "bootstraps             : ${params.bootstrap}"
log.info "experimental design    : ${params.experiment}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

transcriptome_file     = file(params.transcriptome)
exp_file               = file(params.experiment)

/*
 * validate input files
 */
if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing experimental design file: ${exp_file}"

/*
 * Create a channel for read files
 */

Channel
    .fromFilePairs( allReads, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files }


process index {
    input:
    file transcriptome_file

    output:
    file "transcriptome.index" into transcriptome_index

    script:
    //
    // Kallisto tools mapper index
    //
    """
    kallisto index -i transcriptome.index ${transcriptome_file}
    """
}


process mapping {
    tag "reads: $name"

    input:
    file index from transcriptome_index
    set val(name), file(reads) from read_files

    output:
    file "kallisto_${name}" into kallisto_out_dirs, kallisto_out_dirs_viz
    file "stdout.txt" into kallisto_results

    script:
    //
    // Kallisto tools mapper
    //
    def single = reads instanceof Path
    if( !single ) {
        """
        mkdir kallisto_${name}
        kallisto quant -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads} &>stdout.txt
        """
    }
    else {
        """
        mkdir kallisto_${name}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${reads} &>stdout.txt
        """
    }

}


process sleuth {
    input:
    file 'kallisto/*' from kallisto_out_dirs.collect()
    file exp_file

    output:
    file 'sleuth_object.so'
    file 'gene_table_results.txt' into visualisations

    script:
    //
    // Setup sleuth R dependancies and environment
    //

    """
    sleuth.R kallisto ${exp_file}
    """
}


process visualisations {
    publishDir "${params.output}/Visualisations", mode: 'copy'

    container 'lifebitai/vizjson:latest'

    input:
    file mapping from kallisto_out_dirs_viz.collect()
    file gene_table from visualisations

    output:
    file '.report.json' into results

    script:
    """
    tsv2csv.py < $gene_table > gene_table_results.csv
    csv2json.py gene_table_results.csv kallisto 0
    """
}

process multiqc {
    publishDir "${params.output}/MultiQC", mode: 'copy'

    container 'maxulysse/multiqc:1.0'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config
    file ('kallisto/kallisto*') from kallisto_results.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = "--title \"lifebit-ai/kallisto-nf\""
    rfilename = "--filename " + "multiqc_report"
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m kallisto
    """
}
