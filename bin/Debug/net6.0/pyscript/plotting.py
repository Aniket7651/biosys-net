try:
    import sys
    import matplotlib.pyplot as plt
    import seaborn as sns
except ModuleNotFoundError:
    import os
    os.system('pip install -r requirements.txt')

# created a class for reading FASTQ file format
# with takes an argument of file location from the local machine 

################################################ method 2 ###########################################################

def readFASTQfile(fastq_file):
        # seqData is the empty list use for storing data of sequence and ASCII only,
        # seqData variable drop all the sep (+) and header line starting from '@'
        seqData = []
        # seq and ascii list variable stores sequence and thier ASCII code seperately...
        (seq, ascii) = ([], [])
        with open(fastq_file, 'r') as FASTQ:
            # open FASTQ file and read as 'FASTQ', reading start from 2nd lines and store in 'line' variable
            # ie. try to drop first header.
            line = FASTQ.readlines()[1:]
            for i in range(len(line)):
                if i%2 == 0:         # seqData stores those line which contain ASCII and sequence both (even lines)
                    seqData.append(line[i].replace('\n', ''))    # replace all new lines (\n) of the each lines
        
        # this is the seperate part of this program, makes a 'seq' and 'ascii' list from 'seqData', and than 
        # a dict (seqdict) from 'seq' and 'ascii'...
        for i in range(len(seqData)):
            if i%2 == 0:       # this part is same as above even line part, where 'seq' store sequence data from 'seqData'
                seq.append(seqData[i])
            elif i%2 != 0:      # except all (ASCII CODE) stores 'ascii'
                ascii.append(seqData[i])
        seqdict = {key: value for key, value in zip(seq, ascii)}  # use dict. comprehension to make dict of 'seq' and 'ascii'
        # return dict of final sequence and thier respective ASCII code for finding the phred score of that sequence
        return seqdict

def head(fastq_dict, top=5):
        # takes an argument 'top' for reading or returning limited sequence from top 
        # read and select all fastq and makes a list of total tuple items
        lis_limit = list(fastq_dict.items()) 
        return dict(lis_limit[:top])  # and create a dictionary of top limited list

######################################## claculate phred score ########################################

# adding a function to transpose the matrix
def transform(mtrix):
    transf = []  # transposed matrix stored in transf
    for i in range(len(mtrix[0])):  
        lisr = [] 
        for j in range(len(mtrix)):
            lisr.append(mtrix[j][i]) # first store j column of i row in lisr as rows
        transf.append(lisr) # and than lisr store in transf
    return transf

# if our quality score character is in lower case
# so we can use score basded 33 means substract ASCII code to 33
def pred_score_33(str_values):
    # use quality score character as str_values
    pred_list = []
    for char in str_values: 
        # find ascii code and substract with 33 for each char and store these numeric data serise as pred_list
        pred_list.append(ord(char)-33) # ord(character) to find out ASCII code of the perticular char value
    return pred_list

# all the above code is followed for these function also
# but in this case we use based 64 at the place of 33
# based 64 identified by the capital letters only present in the quality score string
def pred_score_64(str_values):
    pred_list = []
    for char in str_values:
        pred_list.append(ord(char)-64) # substract ASCII code with 64 if upper case char only present
    return pred_list

# final function is, to combinly find quality of FASTQ dictionary
def quality(seq_):
    reading_scores = []
    def isBASE():
        base_ = 'ILLUMINA_33'
        for k in seq_:
            ascii = seq_[k]
            for char in ascii:
                if char.islower() == True:
                    base_ = 'ILLUMINA_64'
                break
        return base_

    for k in seq_:
            
        if isBASE() == 'ILLUMINA_33':
            reading_scores.append(pred_score_33(seq_[k]))
        if isBASE() == 'ILLUMINA_64':
                reading_scores.append(pred_score_64(seq_[k]))
    trasf = transform(reading_scores)
    mean = []
    for item in trasf:
        mean.append(sum(item)/len(item))
    return mean, trasf

########################################### FASTQ plotting ############################################

class visualization:

    def GC_graph(self, gc_lis, type='scatter'):
        import matplotlib.pyplot as plt
        import seaborn as sns
    
        Y = [i for i in range(len(gc_lis))]
       
        if type == 'dist' or type == 'hist':
            plt.title('GC CONTENT PER SEQUENCE DISTRIBUTION', font='arial', color='#636057', fontsize=18)
            sns.distplot(gc_lis, kde=True, color='#f78181', bins=30)
            
        else:
            plt.title('GC% PER SEQUENCE READS', font='arial', color='#636057', fontsize=18)
            plt.scatter(gc_lis, Y, color='black', facecolor='#f78181')
            plt.xlabel('GC in %', font='arial', color='#636057', fontsize=12)
            plt.xlim(0, 100)
            plt.ylabel('Number of reads', font='arial', color='#636057', fontsize=12)
        return plt.show()
    
    def scoring_graph_BASE33(self, mean, data, style='default'):

        nt = [i for i in range(len(mean))]
        if style == 'gray':
            # (span1 0-20,  span2 20-28,  span3 28-42,  mean line color,  box color)
            colors = ('#999999', '#cccccc', '#f2f2f2', '#000000', 'white')
        elif style == 'white':
            colors = ('white', 'white', 'white', 'black', '#bfbdbd')
        elif style == 'cool':
            colors = ('#546af7', '#5ca1fa', '#add6f7', '#112ff0', '#f2f207')
        elif style == 'hot':
            colors = ('#f2b750', '#f7cb7e', '#fce4bb', '#f70505', '#f79d8d')
        elif style == 'heatmap':
            colors = ('#65c9f7', '#9ef6f7', '#f79e9e', '#0f0f0f', '#f2a1f7')
        else: 
            colors = ('#ed7272', '#ecfa82', '#82fa8c', '#ff5c33', 'yellow')

        ax = plt.subplots()[1]
        plt.title(f'SCORING GRAPH (ASCII BASED 33), length {len(mean)} bp', size=20, font='arial', color='#636057')
        plt.plot(nt, mean, c=colors[3], linewidth=1)
        sns.boxplot(data, showfliers=False, width=0.9, color=colors[4], linewidth=0.8)
        ranges = [i for i in range(0, len(mean), 5)]
        ax.grid()
        ax.margins(0)
        ax.axhspan(0, 20, facecolor=colors[0], alpha=0.5)
        ax.axhspan(20, 28, facecolor=colors[1], alpha=0.5)
        ax.axhspan(28, 42, facecolor=colors[2], alpha=0.5)
        ax.set_xticks(ranges)
        ax.set_xticklabels(ranges)
        plt.xlim(0, len(nt))
        plt.ylim(0, 42)

        plt.xlabel('Nucleotides (bp)', font='arial', fontsize=12, color='#636057')
        plt.ylabel('Phred Score (Q)', font='arial', fontsize=12, color='#636057')
        plt.savefig('fastq_%s.png'%style)
        return plt.show()

# command lines plotfastq
path = sys.argv[1]
style = sys.argv[2]  # default

fastq_dict = readFASTQfile(path)
# Adict = head(fastq_dict)
qscore = quality(fastq_dict)
visualization().scoring_graph_BASE33(qscore[0], qscore[1], style)

