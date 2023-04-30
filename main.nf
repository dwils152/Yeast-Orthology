#!/users/dwils152/bin nextflow
nextflow.enable.dsl=2

workflow {
    protein_fasta = Channel
                    .fromPath('./data/xyl.faa')
                    .splitFasta(by: 1, file: true )

    //download_db()
    yeast_tax_id()
    blastp(protein_fasta, yeast_tax_id.out)
    to_fasta(blastp.out)
    align_seqs(to_fasta.out)
    build_profile(align_seqs.out)

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
//qseqid = query sequence id
//sseqid = subject sequence id
//sseq = subject sequence (is this only the aligned part?)

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

//process output fasta file from blastp
//since I don't have the data yet, I'll create a dummy file
process align_seqs {
    publishDir "${params.publish_dir}/clustal", mode: 'copy'
    input:
        path fasta
    output:
        path "${fasta.baseName}.aln"
    script:
        """
        clustalo -i ${fasta} -o "${fasta.baseName}.aln"
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
        hmmbuild "${alignment.baseName}.hmm" ${alignment}
        """
}


//process hmm_search {
//
//}

