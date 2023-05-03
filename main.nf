#!/users/dwils152/bin nextflow
nextflow.enable.dsl=2

workflow {
    protein_fasta = Channel
                    .fromPath('./data/xyl.faa')
                    .splitFasta(by: 1, file: true )
    yeast_genomes = Channel
                    .fromPath('./data/yeast_genomes/*.fna')

    //download_db()
    yeast_tax_id()
    blastp(protein_fasta, yeast_tax_id.out)
    to_fasta(blastp.out)
    align_seqs(to_fasta.out)
    build_profile(align_seqs.out)
    get_orfs(yeast_genomes)
    hmm_search(build_profile.out.collect(), get_orfs.out)
    parse_hmm_hits(hmm_search.out)
    hits_per_genomes(parse_hmm_hits.out[0].collect())
    genome_sizes(yeast_genomes.collect())
    plot_size_v_hits(genome_sizes.out, hits_per_genomes.out)
    cat_hits(parse_hmm_hits.out[0].collect())
    //plot_evals(cat_hits.out)
    filter_evals(cat_hits.out)

        yeast_genomes
            .map { [it.name.toString().split("_")[0,1].join(""), it] }
            .set{ genome_map }
            
        parse_hmm_hits.out[0]
            .map { [it.name.toString().split("_")[0,1].join(""), it] }
            .set{ hits_map }
        
        genome_map
            .combine(hits_map, by: 0)
            .map{ id, genome, hits -> [id, genome, hits] }
            .set{ genome_hits }

    extract_fasta(genome_hits)
    cat_fasta(extract_fasta.out.collect())
    compute_kmer_freq(cat_fasta.out)
    align_hmmer_hits(cat_fasta.out)
        
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
        path "${hmmsearch}.tsv", optional: true
        path "${hmmsearch}.bed", optional: true
    script:
        """
        python ${params.scripts}/parse_hmm_hits.py ${hmmsearch}
        """
}

process hits_per_genomes {
    publishDir "${params.publish_dir}/hmmer", mode: 'copy'
    input:
        path all_hits_tables
    output:
        path "hits_per_genome"
    script:
        """
        wc -l * > tmp
        head -n -1 tmp | awk '{print \$1, \$2}' > hits_per_genome
        """
}

process genome_sizes {
    publishDir "${params.publish_dir}/hmmer", mode: 'copy'
    input:
        path genomes
    output:
        path "genome_sizes"
    script:
        """
        genomes=(`echo ${genomes} | tr -d "[],"`)
        for genome in \${genomes[@]}; do
            python ${params.scripts}/fasta_size.py \$genome >> genome_sizes
        done
        """
}

process plot_size_v_hits {
   publishDir "${params.publish_dir}/images", mode: 'copy'
    input:
        path genomes_sizes
        path hmmer_hits
    output:
        path "size_v_hits.png"
    script:
        """
        python ${params.scripts}/plot_size_v_hits.py ${genomes_sizes} ${hmmer_hits}
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
    cache false
    publishDir "${params.publish_dir}/images", mode: 'copy'
    input:
        path hits_table
    output:
        path "e-vals.png"
    script:
        """
        python ${params.scripts}/plot_hmmer_evalues.py ${hits_table}
        """
}

process filter_evals {
    publishDir "${params.publish_dir}/hmmer", mode: 'copy'
    input:
        path hits_table
    output:
        path "${hits_table}.filtered"
    script:
        """
        awk '\$3 <= 0.1' ${hits_table} > ${hits_table}.filtered
        """
}

//add a threshold parameter to the extract_fasta process
process extract_fasta {
    label "Andromeda"
    publishDir "${params.publish_dir}/fasta_hits/single", mode: 'copy'
    input:
        tuple val(id), path(genome), path(hits)
    output:
        path "${id}.faa"
    script:
        """
        python ${params.scripts}/extract_fasta.py ${genome} ${hits} ${id}.faa
        """
}

process cat_fasta {
    publishDir "${params.publish_dir}/fasta_hits", mode: 'copy'
    input:
        path fasta
    output:
        path "all_hits.faa"
    script:
        """
        cat ${fasta} > all_hits.faa
        """
}

process compute_kmer_freq {
    cache false
    publishDir "${params.publish_dir}/kmer_freq", mode: 'copy'
    input:
        path fasta
    output:
        path "${fasta.baseName}.kmer_freq.npy"
    script:
        """
        python ${params.scripts}/kmer_frequencies.py ${fasta} ${fasta.baseName}.kmer_freq.npy
        """
}

process align_hmmer_hits {
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


//align the fasta files
//create a distance matrix
//orthogroup analysis
    //mcl clustering

