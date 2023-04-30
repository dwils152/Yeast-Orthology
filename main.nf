#!/users/dwils152/bin nextflow
nextflow.enable.dsl=2

workflow {
    protein_fasta = Channel
                    .fromPath('./data/xyl.faa')
                    .splitFasta(by: 1, file: true )
    yeast_genomes = Channel
                    .fromPath('./data/**.fsa')

    //download_db()
    yeast_tax_id()
    blastp(protein_fasta, yeast_tax_id.out)
    to_fasta(blastp.out)
    align_seqs(to_fasta.out)
    build_profile(align_seqs.out)
    get_orfs(yeast_genomes)
    hmm_search(build_profile.out.collect(), get_orfs.out)
    parse_hmm_hits(hmm_search.out)
    cat_hits(parse_hmm_hits.out[0].collect())

}

// create a process to download the database... the update_blastdb.pl script comes with the blast+ package but doesn't work
process download_db {
    label "DTN"
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    output:
        path "*"
    script:
        """
        for i in {00..93}; do
            wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.\${i}.tar.gz
            tar -xvzf nr.\${i}.tar.gz
        done
        """
}

// get_species_taxids.sh is a binary that comes with the blast+ package
process yeast_tax_id {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    output:
        path "yeast_tax_id"
    script:
        """
        get_species_taxids.sh -n Saccharomycetes > yeast_tax_id
        """
}

process blastp {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        path protein_fasta
        path yeast_tax_id
    output:
        path "${protein_fasta.baseName}.blastp"
    script:
        """ 
        blastp \
        -query ${protein_fasta} \
        -db nr \
        -out "${protein_fasta.baseName}.blastp" \
        -evalue 10 \
        -outfmt '6 qseqid sseqid sseq' \
        -remote
        """    
}

process to_fasta {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        path blastp_out
    output:
        path "${blastp_out.baseName}.faa"
    script:
        """
        awk 'BEGIN { OFS = "\\n"} {print ">" \$1 \$2, \$3}' ${blastp_out} > "${blastp_out.baseName}.faa"
        """
}

process align_seqs {
    publishDir "${params.publish_dir}/clustal", mode: 'copy'
    input:
        path fasta
    output:
        path "${fasta.baseName}.aln"
    script:
        """
        clustalo -i ${fasta} -o ${fasta.baseName}.aln
        """
}

process build_profile {
    publishDir "${params.publish_dir}/hmmer", mode: 'copy'
    input:
        path alignment
    output:
        path "${alignment.baseName}.hmm"
    script:
        """
        hmmbuild ${alignment.baseName}.hmm ${alignment}
        """
}

process get_orfs {
    publishDir "${params.publish_dir}/orfs", mode: 'copy'
    input:
        path yeast_genomes
    output:
        path "${yeast_genomes.baseName}.orfs.fa"
    script:
        """
        ORFfinder -s 1 -in ${yeast_genomes} -out ${yeast_genomes.baseName}.orfs.fa
        """
}

process hmm_search {
    publishDir "${params.publish_dir}/hmmer/search", mode: 'copy'
    input:
        path hmm
        path orfs
    output:
        path "${orfs.baseName}.hmmsearch"
    script:
        """
        hmmsearch --tblout ${orfs.baseName}.hmmsearch ${hmm} ${orfs} 
        """
}

process parse_hmm_hits {
    publishDir "${params.publish_dir}/hmmer/parsed/tsv", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.publish_dir}/hmmer/parsed/bed", mode: 'copy', pattern: "*.bed"
    input:
        path hmmsearch
    output:
        path "${hmmsearch}.tsv"
        path "${hmmsearch}.bed"
    script:
        """
        python ${params.scripts}/parse_hmm_hits.py ${hmmsearch}
        """
}

process cat_hits {
    publishDir "${params.publish_dir}/hmmer/parsed", mode: 'copy'
    input:
        path hits_table
    output:
        path "all_hits.tsv"
    script:
        """
        cat ${hits_table} > all_hits.tsv
        """
}

process plot_evals {
    publishDir "${params.publish_dir}/images", mode: 'copy'
    input:
        path hits_table
    output:
        path "evals.png"
    script:
        """
        python ${params.scripts}/plot_evals.py ${hits_table}
        """

}

