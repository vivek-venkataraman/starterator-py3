# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Utility functions for Starterator

import pymysql
pymysql.install_as_MySQLdb()
import pymysql as MySQLdb
import configparser as ConfigParser
import getpass
import os, subprocess
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import shutil
# import logging

def decode_if_bytes(value):
    """Handle Python 3 bytes/string conversion for database results"""
    if isinstance(value, bytes):
        return value.decode('utf-8')
    return value


MAKING_FILES = os.path.join(os.path.dirname(os.path.abspath(__file__)))+ "/making_files.py" # absolute path to making files file
ICON_FILE =  os.path.join(os.path.dirname(os.path.abspath(__file__)), "extras", "starterator.svg")
CONFIGURATION_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extras", "starterator.config")
DESKTOP_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extras", "starterator.desktop")
HELP_FILES =  os.path.join(os.path.dirname(os.path.abspath(__file__)), "Help")
PROTEIN_DB =  os.path.join(os.path.dirname(os.path.abspath(__file__)), "Proteins")
INTERMEDIATE_DIR = ""
FINAL_DIR = ""
BLAST_DIR = ""
CLUSTAL_DIR = ""
config_file = os.path.abspath(os.path.join(
    os.getenv("STARTERATOR_CONFIG_DIR", os.path.join(os.environ["HOME"], ".starterator")), 
    "starterator.config"))
icon_file = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator/starterator.svg"))
desktop_file = os.path.abspath(os.path.join(os.environ["HOME"], ".local/share/applications/", "startertor.desktop"))
help_files = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Help")
glade_file = os.path.join(os.path.dirname(os.path.abspath(__file__))) + "/starterator.glade"

class StarteratorError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def get_pham_no(db, phage_name, gene_number):
    """
        Function that gets the phamily number/name of a gene,
        given a geneID from the Mycobacteriophage protein database
        returns a pham number as a string
    """

    # try:
    cursor = db.cursor()
    #Query the database for genes that match the given gene number and phage name
    like_phage_name = phage_name +'%'
    # print like_phage_name, gene_number
    cursor.execute("SELECT `gene`.`GeneID` , `pham`.`Name` , `phage`.`PhageID`\n\
    FROM `gene`\n\
    JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`\n\
    JOIN `phage` ON `gene`.`PhageID` = `phage`.`PhageID`\n\
    WHERE `phage`.`Name` LIKE %s \n\
    AND `gene`.`GeneID` LIKE %s \n\
    ESCAPE '!'", (like_phage_name, '%!_'+ str(gene_number)))
    results = cursor.fetchall()
    # print results, phage_name, gene_number
    # There should only be one result.
    row = results[0]
    pham_no = row[1]
    # print 'pham_no', pham_no
    # Returns the pham number as a string
    return str(pham_no)
    # except:
    #     pass


def find_phams_of_a_phage(db, phage):
    """
        Given a phage in Phamerator, function finds the phams of the genes
        of that phage.
        Returns a list of phams and the length of the nucleotide sequence of the phage 
    """
    # write sql statement to get phams of a phage
    cursor = db.cursor()
    cursor.execute("""SELECT `pham`.`GeneID`, `pham`.`Name`, `phage`.`Name`, `phage`.`SequenceLength`\n\
    from `pham` \n\
    join `gene` on `pham`.`GeneID` = `gene`.`GeneID` \n\
    join `phage` on `phage`.`PhageID` = `gene`.`PhageID`\n\
    where `phage`.`Name` = %s""", (phage))
    results = cursor.fetchall()
    phage_phams = []
    seq_length = results[0][3]
    for row in results:
        # print row[0], row[1]
        phage_phams.append([row[0],str(row[1])])
    return phage_phams, seq_length


def get_protein_sequences():
    proteins = []
    get_db().execute('SELECT GeneID, Translation from gene')
    results = cursor.fetchall()
    for row in results:
        translation = decode_if_bytes(row[1])
        protein = SeqRecord(Seq(translation), id=row[0]+"+", name=row[0], description=row[0])
        proteins.append(protein)
    return proteins

def update_protein_db():
    global BLAST_DIR, PROTEIN_DB
    proteins = get_protein_sequences()
    fasta_file = PROTEIN_DB + 'Proteins'
    count = SeqIO.write(proteins, fasta_file + ".fasta", 'fasta')
    if True:
        blast_db_command = [BLAST_DIR + 'makeblastdb',
                    '-in',"\""+ fasta_file + ".fasta" +"\"",
                    "-dbtype","prot", "-title", "Proteins",
                     "-out", "%s"% fasta_file]
        # print blast_db_command
    # else:
    #     blast_db_command = [BLAST_DIR + 'formatdb',
    #                 '-i', "\""+ fasta_file+ "\"",
    #                 '-o', 'T',
    #                 "-t", "Proteins"]
    #     print blast_db_command
    subprocess.check_call(blast_db_command) #TODO add file clean-up here

def check_protein_db(count):
    get_db().execute('SELECT count(*) from gene')
    results = cursor.fetchall()
    new_count = results[0][0]
    if new_count != count:
        update_protein_db()


def get_gene_number(geneID):
    """
        Get the gene number as an int given a geneID
    """
    # print geneID
    match = re.search(r'^(\w+)([_-]*\w*)_([a-zA-Z]*)(\d+)+$', geneID)
    gene_num = match.groups()[-1]
    # print 'regular exp:', gene_num
    return int(gene_num)

def add_desktop_file():
    '''
    desktop = ConfigParser.RawConfigParser()
    desktop.optionxform = str
    desktop.readfp(open(DESKTOP_FILE))
    desktop_info = dict(desktop.items("Desktop Entry"))
    desktop_info["Exec"] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Starterator")
    desktop_info["Icon"] = os.path.join(os.environ["HOME"], ".starterator/", "starterator.svg")
    desktop.set("Desktop Entry", "Exec",  os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "starterator.sh") )
    desktop.set("Desktop Entry", "Icon", os.path.join(os.environ["HOME"], ".starterator/", "starterator.svg"))
    with open(DESKTOP_FILE, 'wb') as df:
        desktop.write(df)
    subprocess.check_call(["xdg-desktop-icon",  "install", "--novendor", DESKTOP_FILE])
    print("icon file added")
    subprocess.check_call(["xdg-desktop-menu", "install", "--novendor", DESKTOP_FILE])
    print("desktop menu icon installed")
    # if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "starterator.svg")):
    shutil.copyfile(ICON_FILE,
            os.path.join(os.environ["HOME"], ".starterator/", "starterator.svg"))
            '''
         # Function implementation
    pass

def create_folders():
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator"))
    if not os.path.exists(PROTEIN_DB):
        os.mkdir(PROTEIN_DB)
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "starterator.config")):
        if os.path.exists(CONFIGURATION_FILE):
            shutil.copyfile(CONFIGURATION_FILE, 
                os.path.join(os.environ["HOME"], ".starterator", "starterator.config"))
        else:
            # Create default config if template doesn't exist
            config_info = create_default_config()
            write_to_config_file(config_info)
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files"))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Report Files")): 
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Report Files"))

def move_config_file():
    """Create default config file if template doesn't exist"""
    config_path = os.path.join(os.environ['HOME'], '.starterator', 'starterator.config')
    if os.path.exists(CONFIGURATION_FILE):
        shutil.copyfile(CONFIGURATION_FILE, config_path)
        config = get_config()
        config["intermediate_file_dir"] = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files"))
        config["final_file_dir"] = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator", "Report Files"))
        write_to_config_file(config)
    else:
        # Create default config if template doesn't exist
        config_info = create_default_config()
        write_to_config_file(config_info)


def set_up():
    create_folders()
    move_config_file()
    # Comment out or remove the following line to skip setting up desktop files
    # add_desktop_file()

def write_to_config_file(config_info):
    global INTERMEDIATE_DIR, FINAL_DIR, PROTEIN_DB, BLAST_DIR, CLUSTAL_DIR
    config = ConfigParser.RawConfigParser()
    INTERMEDIATE_DIR = config_info["intermediate_file_dir"]
    FINAL_DIR = config_info["final_file_dir"]
    # PROTEIN_DB = config_info["protein_db"]
    CLUSTAL_DIR = config_info["clustalw_dir"]
    BLAST_DIR = config_info["blast_dir"]
    config.add_section('Starterator')
    for name in config_info:
        config.set('Starterator', name, config_info[name])
    print('writing to config file')
    with open(config_file, 'w') as configfile:
        config.write(configfile)

def create_default_config():
    """Create default configuration with environment variable fallbacks"""
    config_dir = os.getenv("CONFIG_DIR", os.path.join(os.environ["HOME"], ".starterator"))
    intermediate_dir = os.getenv("INTERMEDIATE_DIR", os.path.join(config_dir, "Intermediate Files"))
    final_dir = os.getenv("FINAL_DIR", os.path.join(config_dir, "Report Files"))

    return {
        "intermediate_file_dir": intermediate_dir,
        "final_file_dir": final_dir,
        "protein_database": os.path.join(os.path.dirname(os.path.abspath(__file__)), "Proteins"),
        "blast_dir": os.getenv("BLAST_DIR", "/usr/bin/"),
        "clustalw_dir": os.getenv("CLUSTALW_DIR", "/usr/bin/"),
        "database_server": os.getenv("STARTERATOR_DB_HOST", "localhost"),
        "database_name": os.getenv("STARTERATOR_DB_NAME", ""),
        "database_user": os.getenv("STARTERATOR_DB_USER", ""),
        "database_password": os.getenv("STARTERATOR_DB_PASSWORD", ""),
        "database_port": os.getenv("STARTERATOR_DB_PORT", "3306"),
        "count": os.getenv("STARTERATOR_GENE_COUNT", "0")
    }

def get_config():
    """Load user configuration with defaults and environment variable fallbacks"""
    global INTERMEDIATE_DIR, FINAL_DIR, PROTEIN_DB, BLAST_DIR, CLUSTAL_DIR
    config_file = os.path.abspath(os.path.join(
    os.getenv("STARTERATOR_CONFIG_DIR", os.path.join(os.environ.get("HOME", os.environ["HOME"]), ".starterator")), 
    "starterator.config"))
    
    # Use environment variables for configurable paths
    config_dir = os.getenv("STARTERATOR_CONFIG_DIR", os.path.join(os.environ.get("HOME", os.environ["HOME"]), ".starterator"))
    
    # Create directories if they don't exist
    os.makedirs(config_dir, exist_ok=True)
    
    config = ConfigParser.RawConfigParser()
    if not os.path.exists(config_file):
        # Create default configuration
        config_info = create_default_config()
        write_to_config_file(config_info)
    else:
        config.read(config_file)
        config_info = dict(config.items('Starterator'))

    # Use environment variables with config file fallbacks
    INTERMEDIATE_DIR = config_info.get("intermediate_file_dir", os.getenv("STARTERATOR_INTERMEDIATE_DIR", "/tmp/starterator_temp"))
    FINAL_DIR = config_info.get("final_file_dir", os.getenv("STARTERATOR_FINAL_DIR", "/app/reports"))
    CLUSTAL_DIR = config_info.get("clustalw_dir", os.getenv("CLUSTALW_DIR", "/usr/bin/"))
    BLAST_DIR = config_info.get("blast_dir", os.getenv("BLAST_DIR", "/usr/bin/"))
    PROTEIN_DB = os.path.abspath(config_info.get("protein_database", os.path.join(os.path.dirname(os.path.abspath(__file__)), "Proteins")))
    
    # Create needed directories
    os.makedirs(INTERMEDIATE_DIR, exist_ok=True)
    os.makedirs(FINAL_DIR, exist_ok=True)
    
    return config_info

def create_folders():
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator"))
    if not os.path.exists(PROTEIN_DB):
        os.makedirs(PROTEIN_DB, exist_ok=True)
    if not os.path.exists(config_file):
        if os.path.exists(CONFIGURATION_FILE):
            shutil.copyfile(CONFIGURATION_FILE, config_file)
        else:
            # Create default config if template doesn't exist
            config_info = create_default_config()
            write_to_config_file(config_info)
        # logging.info("Created configuration file at {}".format(config_file))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files"))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Report Files")): 
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Report Files"))

def db_connect(config_info):
    db = MySQLdb.connect(config_info['database_server'], 
            config_info['database_user'],
            config_info['database_password'],
            config_info['database_name'])
    return db
    
def attempt_db_connect(config_info):
    try:
        print('attempting to connect', config_info)
        db = MySQLdb.connect(config_info['database_server'], 
                config_info['database_user'],
                config_info['database_password'],
                config_info['database_name'])
        db.close()
    except:
        config_info['database_server'] = input("Enter Database server: ")
        config_info['database_user'] = input("Enter Database username: ")
        config_info['database_password'] = getpass.getpass('Enter Database password: ')
        config_info['database_name'] = input("Enter database name: ")
        print(config_info)
        attempt_db_connect(config_info)
    return config_info

def clean_up_files(file_dir):
    for f in os.listdir(file_dir):
        os.remove(os.path.join(file_dir, f))

