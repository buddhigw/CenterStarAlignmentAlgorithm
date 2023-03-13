#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define GAP_CHAR '-'


int M[200][200];
char Seq[200][200];
int GAP=-5;
int MATCH=10;
int MISSMATCH=-5;
char gaps='-';
char space=' ';

//Reverse the string
char *strrev(char *str){
    char c, *front, *back;

    if(!str || !*str)
        return str;
    for(front=str,back=str+strlen(str)-1;front < back;front++,back--){
        c=*front;*front=*back;*back=c;
    }
    return str;
}

//Step 1: Initialize the table
void initialize(char seq1[],char seq2 []){

   // int seq1len = strlen(seq1);
    //int seq2len = strlen(seq2);
    memset(M, 0, sizeof(M));
}

//Character match score
int match_score(char alpha,char beta)
{
    if(alpha==beta)
    {
        return MATCH;
    }
    else if(alpha==GAP_CHAR || beta==GAP_CHAR)
    {
        return GAP;
    }
    else
    {
        return MISSMATCH;
    }
}

//Step 2: ScoreTable creation
void scoreTable(char seq1[],char seq2 []){
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);

    for(int i=0;i<seq1len+1;i++)
    {
        M[i][0]=GAP*i;
    }

    for(int j=0;j<seq2len+1;j++)
    {
        M[0][j]=GAP*j;
    }

    for(int i=1;i<seq1len+1;i++)
    {
        for(int j=1;j<seq2len+1;j++)
        {
            int match=M[i-1][j-1]+match_score(seq1[i-1],seq2[j-1]);
            int del=M[i-1][j]+GAP;
            int insrt=M[i][j-1]+GAP;
            M[i][j]=MAX(MAX(match,del),insrt);
        }
    }
}

//Getting the pairwise score
int finalize(char seq1[],char seq2[])
{
    int seq1Len = strlen(seq1);
    int seq2Len = strlen(seq2);

    strrev(seq1);
    strrev(seq2);

    int i=0;
    int j=0;
    char symbol[200]="";
    int found=0;
    int score=0;
    float identity=0;

    for(int i=0;i<seq1Len;i++)
    {
        if(seq1[i]==seq2[i])
        {
            strncat(symbol,&seq1[i],1);
            identity+=1;
            score=score+match_score(seq1[i],seq2[i]);
        }
        else if(seq1[i]!=seq2[i] && seq1[i]!=GAP_CHAR && seq2[i]!=GAP_CHAR)
        {
            strncat(symbol,&space,1);
            score=score+match_score(seq1[i],seq2[i]);
            found=0;
        }
        else if(seq1[i]==GAP_CHAR || seq2[i]==GAP_CHAR)
        {
            strncat(symbol,&space,1);
            score=score+GAP;
        }
    }

    identity=identity/seq1Len*100;

    return score;
}

//Get pairwise alignments
struct getAlignments{
        char a1[100];
        char a2[100];
    };
typedef struct getAlignments Struct;

Struct tracebackTable(char seq1[],char seq2 []){
    int i = strlen(seq1);
    int j = strlen(seq2);
    char align1[200]="";
    char align2[200]="";

    while(i>0 && j>0)
    {
        int score_current=M[i][j];
        int score_diagonal=M[i-1][j-1];
        int score_up=M[i][j-1];
        int score_left=M[i-1][j];

        if(score_current==score_diagonal+match_score(seq1[i-1],seq2[j-1]))
        {
            strncat(align1,&seq1[i-1],1);
            strncat(align2,&seq2[j-1],1);
            i--;
            j--;
        }
        else if(score_current==score_left+GAP)
        {
            strncat(align1,&seq1[i-1],1);
            strncat(align2,&gaps,1);
            i--;
        }
        else if(score_current==score_up+GAP)
        {
            strncat(align1,&gaps,1);
            strncat(align2,&seq2[j-1],1);
            j--;
        }
    }

    while(i>0)
    {
        strncat(align1,&seq1[i-1], 1);
        strncat(align2,&gaps,1);
        i--;
    }
    while(j>0)
    {
        strncat(align1,&gaps,1);
        strncat(align2,&seq2[j-1], 1);
        j--;
    }

    strrev(align1);
    strrev(align2);

	Struct s;
    strcpy(s.a1,align1);
    strcpy(s.a2,align2);

    return s;
}

//get string length
int stringLength(char s[])
{
   int c = 0;

   while(s[c] != '\0')
      c++;

   return c;
}

void init_seq()
{
    //char S[6][15] = {"ATCGTGGTACTG","CCGGAGAACTAG","AACGTGCTACTG","ATGGTGAAGTG","CCGGAAAACTTG","TGGCCCTGTATC"};
	//char S[6][8] = {"AGTG","ATCC","ATCG","TCCT","TCGA","TGCG"};
	//char S[3][8] = {"ATCG","ATCC", "ATCG"};
	//char S[2][8] = {"ACAG","ACAAG"};
	char S[5][15] = {"ATTCGGATT","ATCCGGATT","ATGGAATTTT","ATGTTGTT","AGTCAGG"};
	for(int i=0;i<200;i++)
    {
        for(int j=0;j<200;j++)
        {
            Seq[i][j]=S[i][j];
        }
    }

}
int main()
{
	
    init_seq();

    int length = sizeof(Seq)/sizeof(Seq[0]); //get the size of array S

	//create function to this
    int scorematrix[length][length];
    memset(scorematrix, 0, sizeof(scorematrix)); //initialized to 0
    int i, j;
    for(i = 0; i <length; i++)
    {
    	for(j = 0; j <length; j++)
    	{
    		if(i!=j)
    		{
    		    char align1[200]="";
    		    strcat(align1, Seq[i]);
                char align2[200]="";
                strcat(align2, Seq[j]);

    	        initialize(align1,align2);

                scoreTable(align1,align2);

                Struct result;
	            result = tracebackTable(align1,align2);

    			scorematrix[i][j] = finalize(result.a1,result.a2);
    		}
    		else
    		{
    		    scorematrix[i][j] = 0;
    		}
    	}
    }


    int sumOfScore[length];
    int sum=0;
    int centerIndex;

    for(int i=0;i<length;i++)
    {
        sum=0;
        for(int j=0;j<length;j++)
        {
            sum=sum+scorematrix[i][j];
        }
        sumOfScore[i]=sum;
    }



    //find the maximum scored sequence
    int maximum=sumOfScore[0];
    int index=0;
    for(int i=0;i<length;i++)
    {
        if(sumOfScore[i]>maximum)
        {
            maximum=sumOfScore[i];
            index=i;
        }
    }
    printf("Maximum score: %d\n",maximum);
    //find center sequence.
    char *center=Seq[index];
    printf("Center Sequence: %s\n ",center);

    char seqList[length][200];

    //get the alignments with center sequence
    strcpy(seqList[0], center);

    int  indexNum = 1;
    for(int i=0;i<length;i++)
    {
        if(i!=index)
        {
            Struct r;
            initialize(center,Seq[i]);
            scoreTable(center,Seq[i]);
            r=tracebackTable(center,Seq[i]);

            strcpy(seqList[indexNum], r.a2);
            indexNum++;
        }
    }


    char maxLenseq[200];
    int max=0;
    int maxLenSeq = 0;

    //find the mas length sequence
    for(int i=0;i<length;i++)
    {
        if(stringLength(seqList[i])>max)
        {
            strcpy(maxLenseq, seqList[i]);
            max=stringLength(seqList[i]);
            maxLenSeq = i;
        }
    }

    printf("Max len seq: %s\n",maxLenseq);

    char msa[length][200];
    strcpy(msa[0], maxLenseq);

    //Add gaps at the end of the sequences until the length equal to max length sequence.
    int indexJ =1;
	for(int i=0;i<length;i++)
    {
		char strTmp[200];
		strcpy(strTmp, seqList[i]);
        if(i!=maxLenSeq)
        {
			while(max != stringLength(strTmp))
			{
				strncat(strTmp,&gaps,1);
			}
			strcpy(msa[indexJ], strTmp);
            indexJ++;
        }
    }

    //Print all the sequences
    for(j = 0; j <length; j++)
    {
    	printf("%s \n", msa[j]);
    }

    return 0;
}
