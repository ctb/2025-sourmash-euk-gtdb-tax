EUK_K31 = [
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/bilateria-minus-vertebrates.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/eukaryotes-additional.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/eukaryotes-other.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/fungi.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/metazoa-minus-bilateria.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/plants.k31.sig.zip',
    '/group/ctbrowngrp5/2025-genbank-eukaryotes/vertebrates.sig.zip',
    ]

GTDB_K31='/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226-k31.dna.rocksdb'
GTDB_K51='/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226-k51.dna.rocksdb'

DBLIST=[GTDB_K51,
        '/group/ctbrowngrp5/sourmash-db.new/ncbi-euks-2025.01/ncbi-euks-all-2025.01-k51.dna.rocksdb/',
        ]
TAXLIST=['/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.sqldb',
         '/group/ctbrowngrp5/sourmash-db.new/ncbi-euks-2025.01/ncbi-eukaryotes.2025.01.lineages.sqldb',
         ]

SCALED=10_000
KSIZE=31

SAMPLES = ['SRR11125249',
           'SRR12795785',
           'SRR15057925',
           'SRR15057930',
           'SRR5241537']


rule all:
    input:
        expand('{sample}.k31.s10k.gather.with-lineages.csv', sample=SAMPLES),


# do initial gather at k=51
rule gather_sample:
    input:
        metag = "{sample}.trim.sig.zip",
        dblist = DBLIST,
    output:
        gather='{sample}.k51.gather.csv',
        gather_txt='{sample}.k51.gather.out',
    shell: """
        sourmash gather {input.metag} {input.dblist} -k 51 --scaled 10_000 \
           --threshold-bp=0 -o {output.gather} > {output.gather_txt}
    """
    

# then, extract euk matches at k=31
rule collect_k31_euks:
    input:
        gather='{sample}.k51.gather.csv',
        dbs=EUK_K31,
    output:
        '{sample}.euk-k31.matching.sig.zip',
    shell: """
        sourmash sig downsample --picklist {input.gather}:name:ident \
            --scaled {SCALED} -k {KSIZE} \
            {input.dbs} -o {output}
    """
        
        
rule combined_gather_k31_10k:
    input:
        metag = "{sample}.trim.sig.zip",
        gtdb=GTDB_K31,
        euk='{sample}.euk-k31.matching.sig.zip',
    output:
        csv='{sample}.k31.s10k.gather.csv',
        txt='{sample}.k31.s10k.gather.out',
    shell: """
        sourmash gather -k {KSIZE} --scaled {SCALED} --threshold-bp=0 \
            {input.metag} {input.gtdb} {input.euk} \
            -o {output.csv} > {output.txt}
    """

rule annotate_k31_10k:
    input:
        csv='{sample}.k31.s10k.gather.csv',
        tax=TAXLIST,
    output:
        '{sample}.k31.s10k.gather.with-lineages.csv',
    shell: """
        sourmash tax annotate -t {input.tax} -g {input.csv}
    """
