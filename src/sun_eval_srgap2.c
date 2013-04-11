//Xander Nuttle
//sun_eval_srgap2.c
//Evaluates SRGAP2 SUNs with regard to allele and genotype frequencies in individuals having a high-confidence SRGAP2 copy number genotype of 2222.
//
//For each SRGAP2 SUN included in a MIP target region, this program outputs the allele frequency of the SUN (counts of + and - alleles) as well as
//genotype frequencies (counts of +/+, +/-, and -/- genotypes), taking into account only cases where the full genotype can be unambiguously specified.
//Complete genotype specification is necessary for Hardy-Weinberg equilibrium analysis.
//
//Call: ./sun_eval_srgap2 gzipped_mipseqs_file suns_file_with_scores individuals_with_high_confidence_2222_genotypes_file
//Example call: ./sun_eval_srgap2 ~/MIPs/experiments/pos_ctrl_expt/final_results/pos_ctrl_expt.mipseqs.gz ~/MIPs/SRGAP2/SUN_analysis/results/all.suns.scored ~/MIPs/experiments/pos_ctrl_expt/genomes_hc_2222.txt
//
//NOTE: order of genomes listed in high_confidence_genotypes_file must be the order such genomes appear in the gzipped_mipseqs_file!!!
//NOTE: SUNs with scores of -1.000000 have no overlapping SUNKs and thus could not be evaluated using high-coverage genomes (score set to -1.000000 by default)
//
//Strategy:
//	-read in all SUNs so they can be processed repeatedly for different individuals
//	-get 1st set of MIP sequences for 1st individual
//	-get 1st SUN
//	-while master sequence location of 1st SUN < master sequence location of 1st MIP target, get next SUN
//	-if current SUN is in MIP target, analyze it (calculate genotype in that individual)
//	-else get next set of MIP sequences and compare SUN location to location of next MIP target
//
//analyze SUN:
//	-if assignment to a copy is not ambiguous, easy to call SUN as + or - for that copy
//	-if there is some ambiguity in assignment to a copy, all possible assignments to that copy must be SUN + or SUN - for the SUN to be called as + or - in that copy
//	-if both alleles are called, increment appropriate allele counts (+,-) and genotype counts (+/+,+/-,-/-) for that SUN
//	-if at least 1 allele is not called, continue without incrementing any counts
//
//output:
//
//		info to output: paralog SUN_contig_location SUN_master_sequence_location SUN_score Num_SUN_+ Num_SUN_- Num_++ Num_+- Num_--
//		output file: '.suneval'
//

#include<stdio.h>
#include<stdlib.h>
#include<zlib.h>
#include<string.h>
#include<ctype.h>

//set up structure to hold sequence data for a MIP-individual combination
struct mip_data
{
	char indiv[51];
  char gene[51];
  long mip_coord;
  int srgap2_cn_genotype[4];
  char paralog_assigned;
  long contig_mapping_loc;
  char seq[301];
  char qual[301];
  long num_copies_allele;
	long seq_count;
  double prob;
};

//set up structure to hold data for SUNs of interest
struct sun_info
{
	char plog;
	long contig_loc;
	long master_loc;
	double score;
	long num_plus;
	long num_minus;
	long num_pp;
	long num_pm;
	long num_mm;
	int analyzed;
};

//declare functions
void initialize_mip_data(struct mip_data*data);
int compfun(const void*p1,const void*p2);
void collapse(struct mip_data*data,int*nrows);
void reduce(struct mip_data*data,int*nrows);
void analyze_sun(struct sun_info*sundata,struct mip_data*data,int*nrows);
void eval_sun(long sunloc,long seqstart,char paralog,char*sequence,long numcopies,int*p,int*m,int*ap,int*am);
void testsun(char para,char b,long num,int*p,int*m,int*ap,int*am);

int main(int argc,char*argv[])
{
	//set up memory for mip data structure
	struct mip_data*mip_seq_data;
	mip_seq_data=(struct mip_data*)malloc(80*sizeof(struct mip_data)); //allows for lots of lines with ambiguous assignment

	//read in SUNs file with scores and determine the number of SUNs
	FILE*sunsfile;
	sunsfile=fopen(*(argv+2),"r");
	fpos_t pos;
	char contig_name[101];
	long cc,mc;
	double sc;
	long num_suns=0;
	fgetpos(sunsfile,&pos);
	while(fscanf(sunsfile,"%s %ld %ld %lf",contig_name,&cc,&mc,&sc)==4)
	{
		num_suns++;
	}
	fsetpos(sunsfile,&pos);

	//set up memory for sun data structure	
	struct sun_info*sun_data;
	sun_data=(struct sun_info*)malloc(num_suns*sizeof(struct sun_info));

	//read in SUN data and initialze allele and genotype counts to 0
	int sun=0;
	while(fscanf(sunsfile,"%s %ld %ld %lf",contig_name,&(sun_data[sun].contig_loc),&(sun_data[sun].master_loc),&(sun_data[sun].score))==4)
	{
		switch(contig_name[0])
		{
			case '1':	sun_data[sun].plog='A'; break; //1q32.1_final is SRGAP2A
			case 'q':	sun_data[sun].plog='B'; break; //q21_final_oct2011 is SRGAP2B
			case 'p': sun_data[sun].plog='C'; break; //p12_final_oct_2011 is SRGAP2C
			case 'S': sun_data[sun].plog='D'; break; //SRGAP2D_contig is SRGAP2D
		}
		sun_data[sun].num_plus=0;
		sun_data[sun].num_minus=0;
		sun_data[sun].num_pp=0;
		sun_data[sun].num_pm=0;
		sun_data[sun].num_mm=0;
		sun++;
	}
	fclose(sunsfile);

	//read in file having individuals with high confidence 2222 SRGAP2 copy number genotypes
	FILE*indivsfile;
	int g;
	long num_indivs=0,ind=0;
	indivsfile=fopen(*(argv+3),"r");
	fgetpos(indivsfile,&pos);
	while(fscanf(indivsfile,"%s",contig_name)==1)
	{
		num_indivs++;
	}
	fsetpos(indivsfile,&pos);
	char indivs[num_indivs][51];
	for(ind=0;ind<num_indivs;ind++)
	{
		for(g=0;g<51;g++)
		{
			indivs[ind][g]='\0';
		}
	}
	ind=0;
	while(fscanf(indivsfile,"%s",indivs[ind])==1)
  {
    ind++;
  }
	ind=0;
	fclose(indivsfile);

	//set up '.suneval' output file
	FILE*suneval_file;
	char out_name[51];
	int s;
	for(s=0;s<51;s++)
	{
		out_name[s]='\0';
	}
	int num_chars_to_copy;
	char*new_name_start;
	new_name_start=strrchr(*(argv+1),'/')+1;
	if(strrchr(*(argv+1),'/')==NULL)
	{
		new_name_start=*(argv+1);
	}
	num_chars_to_copy=(strchr(new_name_start,'.')-new_name_start);
	if(num_chars_to_copy>42)
	{
		num_chars_to_copy=42;
	}
	strncpy(out_name,new_name_start,num_chars_to_copy);
	strcat(out_name,".suneval");
	suneval_file=fopen(out_name,"w");
	fprintf(suneval_file,"Paralog\tSUN_contig_location\tSUN_master_seq_location\tSUN_score\tNum_SUN_+\tNum_SUN_-\tNum_++\tNum_+-\tNum_--\n");

	//process gzipped mipseqs file
	char ch;
	char*tab_locations1[10],*tab_locations2[10];
	char copystring[5],line1[801],line2[801],current_indiv[51],next_indiv[51];
	int w,x,y,z;
	long current_mip,next_mip;
	long lower=0,upper=0,sunloc=0;
	int do_break=0;
	gzFile*mipseqsfile;
	mipseqsfile=gzopen(*(argv+1),"r");

	while((ch=(gzgetc(mipseqsfile)))!='\n')
		continue; //do nothing with header line
	z=0;
	y=0;
	for(sun=0;sun<num_suns;sun++)
  {
  	sun_data[sun].analyzed=0; //set all suns to the "not yet analyzed" state
  }
	sun=0;
	while((line1[z]=(gzgetc(mipseqsfile)))!='\n')
	{
		if(line1[z]=='\t')
    {
    	tab_locations1[y]=line1+z;
    	y++;
    }
    z++;
  }
  line1[z]='\0';

	while(ch=(gzgetc(mipseqsfile)))
	{
		if(ch==EOF)
		{
			do_break=1;
		}
		x=0;
		z=0;
		if((strcmp(current_indiv,next_indiv)!=0)||(next_mip<current_mip))
		{
			if(strcmp(current_indiv,next_indiv)!=0)
			{
				for(sun=0;sun<num_suns;sun++)
				{
					sun_data[sun].analyzed=0; //for a new individual, reset all suns to the "not yet analyzed" state
				}
			}
			if((ind<num_indivs)&&(strcmp(next_indiv,indivs[ind+1])==0))
			{
				ind++;
			}
			sun=0;
		}
		line2[z]=ch;
		z++;
		current_mip=0;
		next_mip=0;
		strcpy(current_indiv,"indiv");
		strcpy(next_indiv,"indiv");
		initialize_mip_data(mip_seq_data);
		while((current_mip==next_mip)&&(strcmp(current_indiv,next_indiv)==0))
		{
			y=0;
			if(!do_break)
			{
				while((line2[z]=(gzgetc(mipseqsfile)))!='\n')
				{
					if(line2[z]==EOF)
					{
						do_break=1;
						break;
					}
					if(line2[z]=='\t')
      		{
        		tab_locations2[y]=line2+z;
        		y++;
      		}
					z++;
				}
				line2[z]='\0';
			}
			if(mip_seq_data[x].mip_coord=-1)
			{
				strncpy(mip_seq_data[x].indiv,line1,(tab_locations1[0]-line1));
				mip_seq_data[x].indiv[tab_locations1[0]-line1]='\0';
				strncpy(mip_seq_data[x].gene,tab_locations1[0]+1,tab_locations1[1]-tab_locations1[0]-1);
				mip_seq_data[x].gene[tab_locations1[1]-tab_locations1[0]-1]='\0';
				mip_seq_data[x].mip_coord=strtol(tab_locations1[1]+1,NULL,10);
				for(w=0;w<4;w++)
				{
					strncpy(copystring,tab_locations1[2]+1+w,1);
					mip_seq_data[x].srgap2_cn_genotype[w]=atoi(copystring);
				}
				mip_seq_data[x].paralog_assigned=*(tab_locations1[3]+1);
				mip_seq_data[x].contig_mapping_loc=strtol(tab_locations1[4]+1,NULL,10);
				strncpy(mip_seq_data[x].seq,tab_locations1[5]+1,tab_locations1[6]-tab_locations1[5]-1);
				mip_seq_data[x].seq[tab_locations1[6]-tab_locations1[5]-1]='\0';
				strncpy(mip_seq_data[x].qual,tab_locations1[6]+1,tab_locations1[7]-tab_locations1[6]-1);
				mip_seq_data[x].qual[tab_locations1[7]-tab_locations1[6]-1]='\0';
				mip_seq_data[x].num_copies_allele=strtol(tab_locations1[7]+1,NULL,10);
				mip_seq_data[x].seq_count=strtol(tab_locations1[8]+1,NULL,10);
				mip_seq_data[x].prob=strtod(tab_locations1[9]+1,NULL);
				x++;
			}
			if(!do_break)
			{
				current_mip=strtol(tab_locations1[1]+1,NULL,10);
				next_mip=strtol(tab_locations2[1]+1,NULL,10);
				strncpy(current_indiv,line1,(tab_locations1[0]-line1));
				strncpy(next_indiv,line2,(tab_locations2[0]-line2));
				current_indiv[tab_locations1[0]-line1]='\0';
				next_indiv[tab_locations2[0]-line2]='\0';
				strcpy(line1,line2);
				y=0;
				z=0;
				while(line1[z]!='\0')
				{
					if(line1[z]=='\t')
					{
						tab_locations1[y]=line1+z;
						y++;
					}
					z++;
				}
				z=0;
			}
			else
			{
				current_mip=0;
				next_mip=1;
			}
		}

		if(strcmp(current_indiv,indivs[ind])==0)
		{
			qsort(mip_seq_data,80,sizeof(struct mip_data),compfun);
      collapse(mip_seq_data,&x);
      reduce(mip_seq_data,&x);
			lower=current_mip-75;
			upper=current_mip+76;
			if(sun<num_suns)
			{
				sunloc=sun_data[sun].master_loc;
				while(sunloc<lower)
				{
					sun++;
					if(sun==num_suns)
					{
						sunloc=1000000; //a number greater than the master sequence length
						break;
					}
					sunloc=sun_data[sun].master_loc;
				}
			}
			while(sunloc<=upper)
			{
				//printf("Analyzing SUN! indiv=%s SRGAP2%c %ld\n",current_indiv,sun_data[sun].plog,sun_data[sun].master_loc);
				if(!(sun_data[sun].analyzed))
				{
					analyze_sun(&(sun_data[sun]),mip_seq_data,&x);
				}
				sun++;
				if(sun==num_suns)
        {
       		break;
				}
				else
				{	
					sunloc=sun_data[sun].master_loc;
				}
			}
		}
		if(do_break)
		{
			gzclose(mipseqsfile);
			free(mip_seq_data);
			break;
		}
	}
	
	//print output, cleanup, and exit
	for(sun=0;sun<num_suns;sun++)
	{
		fprintf(suneval_file,"SRGAP2%c\t%ld\t%ld\t%lf\t%ld\t%ld\t%ld\t%ld\t%ld\n",sun_data[sun].plog,sun_data[sun].contig_loc,sun_data[sun].master_loc,sun_data[sun].score,sun_data[sun].num_plus,sun_data[sun].num_minus,sun_data[sun].num_pp,sun_data[sun].num_pm,sun_data[sun].num_mm);
	}
	free(sun_data);
	fclose(suneval_file);
	return 0;
}

void initialize_mip_data(struct mip_data*data)
{
	int i,j;
	for(i=0;i<80;i++)
	{
		strcpy(data[i].indiv,"indiv");
		strcpy(data[i].gene,"gene");
		data[i].mip_coord=-1;
		for(j=0;j<4;j++)
		{
			data[i].srgap2_cn_genotype[j]=-1;
		}
		data[i].paralog_assigned='x';
		data[i].contig_mapping_loc=-1;
		strcpy(data[i].seq,"N\0");
		strcpy(data[i].qual,"!\0");
		data[i].num_copies_allele=0;
		data[i].seq_count=0;
		data[i].prob=0.0;
	}
	return;
}

int compfun(const void*p1,const void*p2)
{
	const struct mip_data*a1=p1;
  const struct mip_data*a2=p2;
	if((a1->paralog_assigned)==(a2->paralog_assigned))
	{
		return(strcmp(a1->seq,a2->seq)); //make sure all identical sequences end up in adjacent rows
	}
	else
	{
		return ((a1->paralog_assigned)>(a2->paralog_assigned));
	}
}

void collapse(struct mip_data*data,int*nrows)
{
	int line_num,i;
	int orig_nlines=*nrows;
	for(line_num=1;line_num<orig_nlines;line_num++)
	{
		if((data[line_num].paralog_assigned!='x')&&(data[line_num].paralog_assigned==data[line_num-1].paralog_assigned)&&(strcmp(data[line_num].seq,data[line_num-1].seq)==0))
		{
			data[79]=data[line_num];
			for(i=line_num;i<orig_nlines;i++)
			{
				data[i]=data[i+1];
			}
			(*nrows)--;
			data[line_num-1].num_copies_allele++;
			line_num--;
		}
	}
	return;
}

void reduce(struct mip_data*data,int*nrows)
{
	int i;
	int copiesA,copiesB,copiesC,copiesD;
	copiesA=data[0].srgap2_cn_genotype[0];
	copiesB=data[0].srgap2_cn_genotype[1];
	copiesC=data[0].srgap2_cn_genotype[2];
	copiesD=data[0].srgap2_cn_genotype[3];
	for(i=0;i<*nrows;i++)
	{
		switch(data[i].paralog_assigned)
		{
			case 'A': copiesA-=data[i].num_copies_allele; break;
			case 'a': data[i].num_copies_allele=(data[i].num_copies_allele<copiesA)?data[i].num_copies_allele:copiesA; break;
			case 'B': copiesB-=data[i].num_copies_allele; break;
			case 'b': data[i].num_copies_allele=(data[i].num_copies_allele<copiesB)?data[i].num_copies_allele:copiesB; break;
			case 'C': copiesC-=data[i].num_copies_allele; break;
			case 'c': data[i].num_copies_allele=(data[i].num_copies_allele<copiesC)?data[i].num_copies_allele:copiesC; break;
			case 'D': copiesD-=data[i].num_copies_allele; break;
			case 'd': data[i].num_copies_allele=(data[i].num_copies_allele<copiesD)?data[i].num_copies_allele:copiesD; break;
			default: printf("Invalid paralog!\n"); exit(1);
		}
	}
	return;
}

void analyze_sun(struct sun_info*sundata,struct mip_data*data,int*nrows)
{
	int plus=0,minus=0,ambigplus=0,ambigminus=0;
	int row;
	for(row=0;row<*nrows;row++)
	{
		if(toupper(data[row].paralog_assigned)==sundata->plog)
		{
			eval_sun(sundata->contig_loc,data[row].contig_mapping_loc,data[row].paralog_assigned,data[row].seq,data[row].num_copies_allele,&plus,&minus,&ambigplus,&ambigminus);
		}
	}	
	switch(plus)
	{
		case 2:
						sundata->num_plus+=2;
						sundata->num_pp++;
						break;
		case 1:
						if(minus)
						{
							sundata->num_plus++;
							sundata->num_minus++;
							sundata->num_pm++;
							break;
						}
						else
						{
							if(ambigplus&&ambigminus)
							{
								break;
							}
							else if(ambigplus>ambigminus)
							{
								sundata->num_plus+=2;
								sundata->num_pp++;
								break;
							}
							else
							{
								sundata->num_plus++;
								sundata->num_minus++;
								sundata->num_pm++;
								break;
							}
						}
		case 0:
						if(minus==2)
						{
							sundata->num_minus+=2;
							sundata->num_mm++;
							break;
						}
						else if(minus)
						{
							if(ambigplus&&ambigminus)
              {
                break;
              }
              else if(ambigplus>ambigminus)
              {
                sundata->num_plus++;
								sundata->num_minus++;
                sundata->num_pm++;
                break;
              }
              else
              {
                sundata->num_minus+=2;
                sundata->num_mm++;
                break;
              }
						}
						else
						{
							if(ambigplus&&ambigminus)
							{
								if((ambigplus==1)&&(ambigminus==1))
								{
									sundata->num_plus++;
									sundata->num_minus++;
									sundata->num_pm++;
									break;
								}
								else
								{
									break;
								}
							}
							else if(ambigplus>ambigminus)
							{
								sundata->num_plus+=2;
								sundata->num_pp++;
								break;
							}
							else
							{
								sundata->num_minus+=2;
								sundata->num_mm++;
								break;
							}
						}
	}	//end switch
	sundata->analyzed=1;
	return;
}

void eval_sun(long sunloc,long seqstart,char paralog,char*sequence,long numcopies,int*p,int*m,int*ap,int*am)
{
	char base;	
	while(*sequence)
	{
		switch(*sequence)
		{
			case 's':
								base='-';
								if(seqstart==sunloc)
								{
									testsun(paralog,base,numcopies,p,m,ap,am);
									return;
								}
								seqstart++;
								sequence+=6;
								break;
			case 'i':
								sequence++;
								while(*sequence!='i')
								{
									sequence++;
								}
								sequence++;
								break;
			case 'd':
								sequence++;
								base='-';
								while(*sequence!='d')
								{
									if(seqstart==sunloc)
									{
										testsun(paralog,base,numcopies,p,m,ap,am);
										return;
									}
									sequence++;
									seqstart++;
								}
								sequence++;
			default:
								base='+';
								if(seqstart==sunloc)
								{
									testsun(paralog,base,numcopies,p,m,ap,am);
									return;
								}
								seqstart++;
								sequence++;
								break;
		}
	}	
	return;
}

void testsun(char para,char b,long num,int*p,int*m,int*ap,int*am)
{
	switch(b)
  {
  	case '+':
    					if(isupper(para))
              {
              	*p+=(int)num;
              	break;
             	}
              else
              {
              	*ap+=(int)num;
              	break;
              }
  	case '-':
              if(isupper(para))
              {
              	*m+=(int)num;
                break;
              }
              else
              {
              	*am+=(int)num;
              }
  }
	return;
}
