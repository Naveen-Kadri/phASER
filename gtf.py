"""
"""

class gtf () :
    infos_to_write=["gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_biotype",   "exon_number", "strand",  "phase",   "start", "end"]
    def __init__(self, line:str) -> None:
        seqid, source, ftype, start,end, score,strand,phase,attrs = line.rstrip().split ("\t")
        self.seqid = seqid
        self.source = source
        self.ftype=ftype
        self.start=int (start)
        self.end=int(end)
        self.score = score
        self.strand=strand
        self.phase=phase
        self.attrs = attrs

        
        for el in attrs.split ("; "):
            spl=el.split(" ")
            if len (spl) ==2:
                atr, val=spl
                setattr(self, atr, val.replace('"', ''))
            else:
                setattr(self, atr, "NA")

    def __str__ (self) ->str :
        return (f"{self.ftype} for the gene {self.gene_id} ({self.gene_biotype})  {self.strand} strand of {self.seqid} from {self.start} to {self.end} ")

                
    @property
    def length(self):
        return abs(self.end - self.start) +1

    @classmethod
    def set_infos_to_write(cls, mylist=infos_to_write):
        cls.infos_to_write=mylist
        
    def print_info(self):
        tw=[]
        for info in gtf.infos_to_write:
            try:
                tw.append (getattr (self, info)   )
            except(AttributeError):
                tw.append ("NA" )

        return (tw)
        
    #overloding -- special/magic function
    def __len__(self):
        return abs(self.end - self.start) +1

    @staticmethod
    def get_introns(exon_list, strand):
        nex=len(exon_list) -1
        intron_list=[]
        if strand == "+":
            for i in range (0, nex):
                intron_list.append ([   exon_list [i] [1] +1, exon_list [i+1] [0] - 1   ])
        else:
            for i in range (0, nex):
                intron_list.append ([   exon_list [i] [1] -1, exon_list [i+1] [0] + 1   ])
        return (intron_list)

    @staticmethod
    def count_var (start, end, phase,varlist):
        if start > end:
            start, end = end,start
        nvar=dict ()
        nbases=dict ()
        tracker=int (phase) -1
        for pos in range (start, end+1):
            tracker+=1
            triplet_pos=tracker%3
            nbases [triplet_pos]= nbases.get(triplet_pos,0)+1
            if pos in varlist:
                nvar [triplet_pos]= nvar.get(triplet_pos,0)+1
        tr=[]
        for k in range (3):
            # 0 is for 3rd pos of the triplet
            tr.append (str  (nbases.get(k,0)   )   )
            tr.append (str   (nvar.get (k,0))   )
        return ("\t".join (tr))
            










