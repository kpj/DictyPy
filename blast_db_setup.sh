set -e
set -u

PATH="$PATH:/home/kpj/blast/ncbi-blast-2.2.30+-src/c++/ReleaseMT/bin"
DB_DIR="/home/kpj/blast/db/nr/"


function download_db {
    (
        cd "$DB_DIR"
        db="${1:-nr}"

        if [[ -e "$db.00.psq" ]] ; then
            echo "Looks like archive was already downloaded, skipping..."
            exit
        fi

        update_blastdb.pl --decompress "$db"
    )
}

function extract_taxid {
    #It might be necessary to create extra swap space:
    #    $ fallocate -l 20G swapfile
    #    $ mkswap swapfile
    #    $ sudo swapon swapfile
    #    ..
    #    $ sudo swapoff -a
    #    $ rm -f swapfile

    taxid="$1"
    append_fname="$2"
    db="${3:-nr}"

    fname="results/blast_db_overview.txt"
    if [[ ! -e "$fname" ]] ; then
        blastdbcmd \
            -db "$DB_DIR/$db" \
            -entry all \
            -outfmt "%g %T" \
            -out "$fname"
    fi

    grep " $taxid$" "$fname" | awk '{print $1}' | tee -a "$append_fname"
}

function create_new_db {
    new_db_name="$1"
    gi_file="$2"
    db="${3:-nr}"

    tmp_file="tmp.fa"

    # extract needed sequences as fasta
    blastdbcmd \
        -db "$DB_DIR/$db" \
        -entry_batch "$gi_file" \
        -out "$tmp_file"

    # create blast database
    makeblastdb \
        -dbtype prot \
        -in "$tmp_file" \
        -input_type fasta \
        -out "$new_db_name" \
        -title "$new_db_name"
}


new_db_name="custom_blast_db"
gi_list="gi_list.txt"

download_db

rm -f "$gi_list"
while read tid ; do
    extract_taxid "$tid" "$gi_list"
done < "taxids.txt"

create_new_db "$new_db_name" "$gi_list"


<<'COMMENT'
# list dbs in current directory
$ blastdbcmd -list . -recursive -list_outfmt '%t (%p): %n seqs'

# export db as fasta file
$ blastdbcmd -db "$db" -entry all -outfmt '%f' -out "$db.fa"

# get info about db
$ blastdbcmd -db "$db" -info
COMMENT
