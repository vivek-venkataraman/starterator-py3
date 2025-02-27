from .database import DB, get_db
from . import phamgene
# don't want mutliples of phage object
# before making a phage, 
phage_list = {}
def new_phage(phage_id=None, name=None, cluster=None, sequence=None):
    if not phage_id:
        phage = Phage(phage_id, name, cluster, sequence)
        phage_id = phage.get_id()
    if phage_list.get(phage_id, True):
        phage_list[phage_id] = Phage(phage_id, name, cluster, sequence)
    return phage_list[phage_id]

class Phage(object):
    def __init__(self, phage_id=None, name=None, cluster=None, sequence=None, phamerated=True):
        self.phamerated = phamerated
        self.name = name
        self.phage_id = phage_id
        self.sequence = sequence
        self.cluster= cluster
        self.genes = None
        self.phams = None
        self.genes = None

    def get_name(self):
        if not self.name:
            row = get_db().get(
                "SELECT Name, Cluster, Sequence from phage where PhageID = %s",
                self.phage_id)
            self.name = row[0]
            self.cluster = row[1]
            self.sequence = row[2]
        return self.name

    def get_id(self):
        if not self.phage_id:
            row = get_db().get(
                "SELECT PhageID, Cluster, Sequence from phage where Name like %s",
                self.name)
            self.phage_id = row[0]
            self.cluster = row[1]
            self.sequence = row[2]
        return self.phage_id
    
    def get_sequence(self):
        if not self.sequence:
            if self.phage_id:
                row = get_db().get(
                    "SELECT Sequence from phage where phageID = %s", self.phage_id)
                self.sequence = row[0]
            elif self.name:
                row = get_db().get(
                    "SELECT Sequence from phage where Name like %s", self.name)
                self.sequence = row[0]
        return self.sequence

    def length(self):
        if not self.sequence:
            self.get_sequence()
        return len(self.sequence)

    def get_cluster(self):
        if not self.cluster:
            if self.phage_id:
                row = get_db().get(
                    "SELECT Cluster from phage where PhageID = %s", self.phage_id)
                self.cluster = row[0]
            elif self.name:
                row = self.db.get(
                    "SELECT Cluster from phage where Name like %s", self.name)
                self.cluster = row[0]
        return self.cluster

    def get_genes(self):
        if not self.genes:
            if not self.phage_id:
                self.get_phage_id()
            self.genes = []
            results = get_db().query(
                "SELECT `pham`.`GeneID`, `pham`.`name`, `gene`.Name, \n\
                `gene`.`Start`, `gene`.`Stop`, `gene`.`Orientation`\n\
                FROM `pham` JOIN `gene` on `pham`.`GeneID` = `gene`.`GeneID`\n\
                WHERE `gene`.`PhageID` = %s", self.phage_id) 
            for row in results:
                gene = phamgene.PhamGene(row[0], row[3], row[4], row[5], pham_no=row[2])
                self.genes.append(gene)
        return self.genes


    def get_phams(self):
        if not self.phams:
            self.get_name()
            self.phams = {}
            # gene.Name can be in from gp<Number>, gene<Number>, or <PHAGE_NAME>_<Number>
            results = get_db().query(
                "SELECT `pham`.`GeneID`, `pham`.`name`, `gene`.Name,\n\
                `gene`.`Start`, `gene`.`Stop`, `gene`.`Orientation`\n\
                FROM `pham` JOIN `gene` on `pham`.`GeneID` = `gene`.`GeneID`\n\
                WHERE `gene`.`PhageID` = %s", self.phage_id)
            for row in results:
                if row[1] not in self.phams:
                    self.phams[row[1]] = []
                gene = phamgene.PhamGene(row[0], row[3], row[4], row[5], self.phage_id)
                self.phams[row[1]].append(gene)
        return self.phams

class UnPhamPhage(Phage):
    def __init__(name, fasta_file, profile_file):
        pass




