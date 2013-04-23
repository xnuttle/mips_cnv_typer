//Xander Nuttle
//mipmaker_score.c
//Selects regions from a multiple alignment of paralogous sequences desirable for MIP targeting for paralog-specific copy number genotyping, also outputting scores
//This specific program is only applicable for SRGAP2 but could be easily modified for other loci
//Call: ./mipmaker_score SRGAP2A_shared.fasta SRGAP2B_shared.fasta SRGAP2C_shared.fasta SRGAP2D_shared.fasta regions_to_target.bed
//
//Scoring system: +100 for 1st SUN in a paralog's 112-mer, +20 for 2nd SUN in a paralog's 112-mer, +5 for 3rd SUN in a paralog's 112-mer, +1 for additional SUNs in a paralog's 112-mer


//GENERALIZATIONS TO MAKE
//-allow for a variable number of aligned paralog sequences to be input
//-allow for aligned paralog sequences > 300 kbp
//-generate variable number of output files, one for MIP targets distinguishing 1 paralog, one for MIP targets distinguinshing 2 paralogs, ..., one for MIP targets distinguishing N-2 paralogs, and one for MIP targets distinguishing N paralogs
//-eliminate scoring based on SUNs (separating output files is better than SUN scoring)
//-eliminate SRGAP2-specific features (such as code related to SRGAP2D internal deletion)


#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc,char*argv[])
{
	FILE*in1,*in2,*in3,*in4,*out;
	long seqsize=300000,i=0,j;
	char*A,*B,*C,*D; //SRGAP2A,SRGAP2B,SRGAP2C,SRGAP2D
	A=(char*)malloc((seqsize+1)*sizeof(char));
	B=(char*)malloc((seqsize+1)*sizeof(char));
	C=(char*)malloc((seqsize+1)*sizeof(char));
	D=(char*)malloc((seqsize+1)*sizeof(char));
	if((A==NULL)||(B==NULL)||(C==NULL)||(D==NULL))
	{
		printf("Memory allocation failed!\n");
		return 1;
	}
	
	//read sequences into character arrays
	in1=fopen(*(argv+1),"r");
	in2=fopen(*(argv+2),"r");
	in3=fopen(*(argv+3),"r");
	in4=fopen(*(argv+4),"r");
	char base;

	while((base=getc(in1))!='\n')
		continue;
	while((base=getc(in2))!='\n')
		  continue;
	while((base=getc(in3))!='\n')
		  continue;
	while((base=getc(in4))!='\n')
		  continue;	
	
	while((base=getc(in1))!=EOF)
	{
		if(!isspace(base))
		{
			A[i]=base;
			i++;
		}
	}	
	
	i=0;
	while((base=getc(in2))!=EOF)
	{
		if(!isspace(base))
		{
			B[i]=base;
			i++;
		}
	}

	i=0;
	while((base=getc(in3))!=EOF)
	{
		if(!isspace(base))
		{
			C[i]=base;
			i++;
		}
	}

	i=0;
	while((base=getc(in4))!=EOF)
	{
		if(!isspace(base))
		{
			D[i]=base;
			i++;
		}
	}
	
	A[i]='\0';
	B[i]='\0';
	C[i]='\0';
	D[i]='\0';
	
	//scan character arrays for suitable target regions
	out=fopen(*(argv+5),"w");
	int target_size=112,arm_length=20,x,sunsA=0,sunsB=0,sunsC=0,sunsD=0;
	long delbases=0,score=0,SRGAP2D_delstart=105947,SRGAP2D_delend=213356;
	char s1[target_size+1],s2[target_size+1],s3[target_size+1],s4[target_size+1];
	s1[target_size]='\0';
	s2[target_size]='\0';
	s3[target_size]='\0';
	s4[target_size]='\0';
	for(i=0;i<=strlen(A)-(target_size+(2*arm_length));i++)
	{
		sunsA=0;
		sunsB=0;
		sunsC=0;
		sunsD=0;
		score=0;
		if(A[i]=='-')
		{
			delbases++;
		}
		if(((i-delbases+1)>=SRGAP2D_delstart)&&((i-delbases+1)<=SRGAP2D_delend))
		{
			if((strncmp(A+i,B+i,arm_length)==0)&&(strncmp(A+i,C+i,arm_length)==0)&&(strncmp(B+i,C+i,arm_length)==0))
			{
				if((strncmp(A+arm_length+target_size+i,B+arm_length+target_size+i,arm_length)==0)&&(strncmp(A+arm_length+target_size+i,C+arm_length+target_size+i,arm_length)==0)&&(strncmp(B+arm_length+target_size+i,C+arm_length+target_size+i,arm_length)==0))
				{
					if((strncmp(A+arm_length+i,B+arm_length+i,target_size)!=0)&&(strncmp(A+arm_length+i,C+arm_length+i,target_size)!=0)&&(strncmp(B+arm_length+i,C+arm_length+i,target_size)!=0))
					{
				
						strncpy(s1,A+i+arm_length,target_size);
						strncpy(s2,B+i+arm_length,target_size);
						strncpy(s3,C+i+arm_length,target_size);
						for(x=0;x<target_size;x++)
						{
							if((strchr(s1,'-')==NULL)&&(strchr(s2,'-')==NULL)&&(strchr(s3,'-')==NULL))
							{
								for(x=0;x<target_size;x++)
            		{
              		if((s1[x]!=s2[x])&&(s1[x]!=s3[x]))
                		sunsA++;
              		if((s2[x]!=s1[x])&&(s2[x]!=s3[x]))
                		sunsB++;
              		if((s3[x]!=s1[x])&&(s3[x]!=s2[x]))
                		sunsC++;
            		}
            		while(sunsA!=0)
            		{
              		if(sunsA==1)
                		score+=100;
              		else if(sunsA==2)
                		score+=20;
              		else if(sunsA==3)
                		score+=5;
              		else
                		score+=1;
              		sunsA--;
            		}
            		while(sunsB!=0)
            		{
              		if(sunsB==1)
                		score+=100;
              		else if(sunsB==2)
                		score+=20;
              		else if(sunsB==3)
                		score+=5;
              		else
                		score+=1;
              		sunsB--;
            		}
            		while(sunsC!=0)
            		{
              		if(sunsC==1)
                		score+=100;
              		else if(sunsC==2)
                		score+=20;
              		else if(sunsC==3)
                		score+=5;
              		else
                		score+=1;
              		sunsC--;
            		}
								fprintf(out,"SRGAP2_1q32\t%ld\t%ld\t%ld\n",i+arm_length-delbases,i+arm_length+target_size-delbases,score);
								//printf("%s\n%s\n%s\n%s\n",s1,s2,s3,s4);
							}			
						}
					}
				}
			}
		}
		else
		{
			if((strncmp(A+i,B+i,arm_length)==0)&&(strncmp(A+i,C+i,arm_length)==0)&&(strncmp(A+i,D+i,arm_length)==0)&&(strncmp(B+i,C+i,arm_length)==0)&&(strncmp(B+i,D+i,arm_length)==0)&&(strncmp(C+i,D+i,arm_length)==0))
			{
				if((strncmp(A+arm_length+target_size+i,B+arm_length+target_size+i,arm_length)==0)&&(strncmp(A+arm_length+target_size+i,C+arm_length+target_size+i,arm_length)==0)&&(strncmp(A+arm_length+target_size+i,D+arm_length+target_size+i,arm_length)==0)&&(strncmp(B+arm_length+target_size+i,C+arm_length+target_size+i,arm_length)==0)&&(strncmp(B+arm_length+target_size+i,D+arm_length+target_size+i,arm_length)==0)&&(strncmp(C+arm_length+target_size+i,D+arm_length+target_size+i,arm_length)==0))
				{
					if((strncmp(A+arm_length+i,B+arm_length+i,target_size)!=0)&&(strncmp(A+arm_length+i,C+arm_length+i,target_size)!=0)&&(strncmp(A+arm_length+i,D+arm_length+i,target_size)!=0)&&(strncmp(B+arm_length+i,C+arm_length+i,target_size)!=0)&&(strncmp(B+arm_length+i,D+arm_length+i,target_size)!=0)&&(strncmp(C+arm_length+i,D+arm_length+i,target_size)!=0))
					{
			     	strncpy(s1,A+i+arm_length,target_size);
						strncpy(s2,B+i+arm_length,target_size);
						strncpy(s3,C+i+arm_length,target_size);
						strncpy(s4,D+i+arm_length,target_size);
						if((strchr(s1,'-')==NULL)&&(strchr(s2,'-')==NULL)&&(strchr(s3,'-')==NULL)&&(strchr(s4,'-')==NULL))
						{
							for(x=0;x<target_size;x++)
              {
              	if((s1[x]!=s2[x])&&(s1[x]!=s3[x])&&(s1[x]!=s4[x]))
                	sunsA++;
                if((s2[x]!=s1[x])&&(s2[x]!=s3[x])&&(s2[x]!=s4[x]))
                  sunsB++;
                if((s3[x]!=s1[x])&&(s3[x]!=s2[x])&&(s3[x]!=s4[x]))
                  sunsC++;
								if((s4[x]!=s1[x])&&(s4[x]!=s2[x])&&(s4[x]!=s3[x]))
									sunsD++;
              }
              while(sunsA!=0)
              {
              	if(sunsA==1)
                	score+=100;
                else if(sunsA==2)
                  score+=20;
                else if(sunsA==3)
                  score+=5;
                else
                  score+=1;
                sunsA--;
              }
              while(sunsB!=0)
              {
              	if(sunsB==1)
                	score+=100;
                else if(sunsB==2)
                  score+=20;
                else if(sunsB==3)
                  score+=5;
                else
                  score+=1;
                sunsB--;
              }
              while(sunsC!=0)
              {
                if(sunsC==1)
                  score+=100;
                else if(sunsC==2)
                  score+=20;
                else if(sunsC==3)
                  score+=5;
                else
                  score+=1;
                sunsC--;
              }
							while(sunsD!=0)
              {
                if(sunsD==1)
                  score+=100;
                else if(sunsD==2)
                  score+=20;
                else if(sunsD==3)
                  score+=5;
                else
                  score+=1;
                sunsD--;
              }
							fprintf(out,"SRGAP2_1q32\t%ld\t%ld\t%ld\n",i+arm_length-delbases,i+arm_length+target_size-delbases,score);
							//printf("%s\n%s\n%s\n%s\n",s1,s2,s3,s4);
						}
					}
				}
			}
		}
	}
	
	//clean up and exit
	fclose(in1);
	fclose(in2);
	fclose(in3);
	fclose(in4);
	fclose(out);
	free(A);
	free(B);
	free(C);
	free(D);
	return 0;
}
