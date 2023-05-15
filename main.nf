#!/users/dwils152/bin nextflow
nextflow.enable.dsl=2

workflow {
    protein_fasta = Channel
                    .fromPath('./data/xyl.faa')
                    .splitFasta(by: 1, file: true )
    yeast_genomes = Channel
                    .fromPath('./data/yeast_genomes/*.fna')
    ko_annotations = Channel
                    .fromPath('./data/ko_annotations.txt')

    //download_db()
    //yeast_tax_id()
    //blastp(protein_fasta)
    //get_subtree_id(yeast_tax_id.out)
    filter_by_taxid(blastp.out, get_subtree_id.out)
    sample_fasta(get_subtree_id.out, filter_by_taxid.out)
    to_fasta(sample_fasta.out)
    align_seqs(to_fasta.out)
    build_profile(align_seqs.out)
    query_length(protein_fasta)
    get_orfs(yeast_genomes, query_length.out.first()) //first converts queue to value
    hmm_search(build_profile.out.collect(), get_orfs.out)
    parse_hmm_hits(hmm_search.out)
    hits_per_genomes(parse_hmm_hits.out[0].collect())
    genome_sizes(yeast_genomes.collect())
    plot_size_v_hits(genome_sizes.out, hits_per_genomes.out)
    cat_hits(parse_hmm_hits.out[0].collect())
    plot_evals(cat_hits.out)
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

    extract_fasta(genome_hits).set { fna }
    fna_to_faa(extract_fasta.out).set { faa }
    cat_fasta(fna.collect(), faa.collect())
    compute_kmer_freq(cat_fasta.out[0])
    cluster_seqs(compute_kmer_freq.out)
    get_hits_w_annot(ko_annotations, cat_fasta.out[1])
    plot_seq_length(cat_fasta.out[1], get_hits_w_annot.out)
    //align_hmmer_hits(get_hits_w_annot.out)
    //run_orthofinder(get_hits_w_annot.out.splitFasta(by: 1, file: true).collect())
    //generate_ml_tree(align_hmmer_hits.out)
    */
        
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
        env id
    script:
        """
        get_species_taxids.sh -n Saccharomycetales > yeast_tax_id
        id=(`head -n2 yeast_tax_id | tail -n1`)
        id=\${id[2]}
        """
}

process blastp {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        path protein_fasta
    output:
        path "${protein_fasta.baseName}.blastp"
    script:
        """ 
        blastp \
        -query ${protein_fasta} \
        -db nr \
        -out "${protein_fasta.baseName}.blastp" \
        -outfmt "6 qseqid sseqid sseq staxid evalue" \
        -remote
        """    
}

process get_subtree_id {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        val yeast_tax_id
    output:
        path "yeast_subtree_id"
    script:
        """
        taxonkit list -i ${yeast_tax_id} | taxonkit reformat -I 1 > yeast_subtree_id
        """
}

process filter_by_taxid {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        path blastp_out
        val yeast_tax_id
    output:
        path "${blastp_out}.filtered"
    script:
        """
        python ${params.scripts}/filter_blast_by_tax.py ${blastp_out} ${yeast_tax_id} "${blastp_out}.filtered"
        """
}

// choose a maximum of 5 sequences per family to build profile
process sample_fasta {
    publishDir "${params.publish_dir}/blastp", mode: 'copy'
    input:
        path yeast_tax
        path fasta
    output:
        path "${fasta}.sampled"
    script:
        """
        python ${params.scripts}/sample_blast_hits.py ${fasta} ${yeast_tax} "${fasta}.sampled"
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

process query_length {
    input:
        path protein_fasta
    output:
        env length
    script:
        """
        length=`grep -v ">" ${protein_fasta} | tr -d '\\n' | wc -c`
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
        ORFfinder -ml \$threshold -in ${yeast_genomes} -out ${yeast_genomes.baseName}.orfs.fa
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

process filter_hits_by_len {
    publishDir "${params.publish_dir}/hmmer/parsed", mode: 'copy'
    input:
        path hits_table
    output:
        path "${hits_table}.filtered"
    script:
        """
        
        """
}

process plot_evals {
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
    publishDir "${params.publish_dir}/fasta_hits/single_fna", mode: 'copy'
    input:
        tuple val(id), path(genome), path(hits)
    output:
        path "${id}.fna"
    script:
        """
        python ${params.scripts}/extract_fasta.py ${genome} ${hits} ${id}.fna
        """
}

process fna_to_faa {
    publishDir "${params.publish_dir}/fasta_hits/single_faa", mode: 'copy'
    input:
        path fna
    output:
        path "${fna.baseName}.faa"
    script:
        """
        python ${params.scripts}/fna_to_faa.py ${fna} ${fna.baseName}.faa
        """
}

process cat_fasta {
    publishDir "${params.publish_dir}/fasta_hits", mode: 'copy'
    input:
        path fna
        path faa
    output:
        path "all_hits.fna"
        path "all_hits.faa"
    script:
        """
        cat ${fna} > all_hits.fna
        cat ${faa} > all_hits.faa
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

process cluster_seqs {
    label "big_mem"
    publishDir "${params.publish_dir}/clustering", mode: 'copy'
    input:
        path fasta
    output:
        path "${fasta.baseName}.tsne.png"
        path "${fasta.baseName}.knn.png"
    script:
        """
        python ${params.scripts}/knn.py ${fasta}
        """
}

process plot_seq_length {
    publishDir "${params.publish_dir}/images", mode: 'copy'
    input:
        path all_hits
        path hits_w_annot
    output:
        path "seq_lengths.png"
    script:
        """
        python ${params.scripts}/plot_seq_lengths.py ${all_hits} ${hits_w_annot} seq_lengths.png
        """
}

process get_hits_w_annot {
    publishDir "${params.publish_dir}", mode: 'copy'
    input:
        path ko_annotations
        path all_hits
    output:
        path "all_hits_w_annot.faa"
    script:
        """
        awk 'NF==2' ${ko_annotations} > w_annot
        python ${params.scripts}/get_hits_w_annot.py ${all_hits} w_annot all_hits_w_annot.faa
        """
}

process align_hmmer_hits {
    label "big_mem"
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

process run_orthofinder {
    publishDir "${params.publish_dir}/orthofinder", mode: 'copy'
    input:
        path fasta
    output:
        "*"
    script:
        """
        mkdir -p run_dir
        mv ${fasta} run_dir
        python ~/OrthoFinder_source/orthofinder.py -f run_dir
        """
}

process generate_ml_tree {
    publishDir "${params.publish_dir}/trees", mode: 'copy'
    input:
        path alignment
    output:
        path "${alignment}.tree"
    script:
        """
        FastTreeDbl ${alignment} > ${alignment}.tree
        """
}



