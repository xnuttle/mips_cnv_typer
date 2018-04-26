//Xander Nuttle
//call_mip_pscn.c
//Call: ./call_mip_pscn miptargets_file mipcounts_file (long)max_pscn
//
//The miptargets file input must be in format v3 (where MIPs have letter-based specificities).

#include<stdio.h>
#include<string.h>
#include<gsl/gsl_randist.h>
#include<float.h>
#include<math.h>
#include<stdlib.h>
#define KRED "\x1B[31m"
#define KYEL "\x1B[33m"
#define L_PRIOR -7.5 //log likelihood increment for each copy difference from copy number 2 (used in calculating prior probabilities for each possible copy number state)
#define L_MIN -30.0 //minimum log likelihood value for any single MIP probe arbitrarily assigned to -30.0
#define L_LRT 40.0 //minimum difference in log-likelihoods between a path with more copy number states and a path with fewer copy number states for caller to consider the former path
#define M_MIN 5 //minimum number of MIPs in each copy number state for caller to consider a path having multiple copy number states
#define LOD_MAX 1000.0 //maximum value of LOD score used to quantify how much more likely the maximally likely 0-transition path is than the next most likely such path

//set up node structure for dynamic programming to find maximum likelihood path through graph
//allow detection of two copy number transitions across spatial extent of targeted sequence (e.g., to detect an internal duplication, deletion, or interlocus gene conversion signature)
struct node
{
	double likelihood;
	double max_0t;
	double max_1t;
	double max_2t;
	struct node*path_1t;
	struct node*path_2t;
};

void count_targets(FILE*miptargs,long*nummips,long*numseqs);
void get_mip_info(FILE*miptargs,long nummips,long speclength,long targlocs[nummips],char specvecs[nummips][speclength]);
void set_priors(long numseqs,long parastates,long numstates,double priorvec[numstates],long cstates[numstates][numseqs]);
void init_output(FILE*outfiles[3],char*base,long numseqs);
void init_graph(struct node*graph,long nummips,long numstates);
void fill_graph(struct node*graph,FILE*mcounts,long numseqs,long nummips,long numstates,long speclength,char specvecs[nummips][speclength],long cnstates[numstates][numseqs]);
void init_counts(long nseqs,unsigned int rcounts[nseqs],unsigned int fcounts[nseqs]);
void init_probs(long nseqs,double pvec[nseqs]);
long num_copies(long state,long nseqs,long nstates,long copystates[nstates][nseqs]);
void parse_graph(struct node*graph,long numseqs,long nummips,long numstates,double priorvec[numstates],long cnstates[numstates][numseqs]);
double segmax(struct node*lgraph,long index,long nseqs,long nstates,long copystates[nstates][nseqs],int path);
int trans_good(long index,long index2,long n_states,long n_seqs,long pscn_states[n_states][n_seqs]);
void get_path_maxes(struct node*lgraph,long nmips,long nstates,double*m0,double*nm0,double*m1,double*m2,long*im0,long*im1,long*im2);
int diff_support(struct node*lgraph,long index,long index2,long nstates);
int assess_path(struct node*lgraph,double m0,double m1,double m2,long im1,long im2,long nmips,long nstates,long*cnstates,long*edgemips,int path);
int imax(int i1,int i2);
double dmax(double d1,double d2);
void print_output(FILE*outfiles[3],char*individual,double m0,double nm0,double m1,double m2,long im0,int ntrans,long*cnstates,long*edgemips,long nummips,long targlocs[nummips],long nstates,long nseqs,long copystates[nstates][nseqs]);
double dmin(double d1,double d2);
void print_pscn(FILE*out,long state,long n_states,long n_seqs,long pscn_states[n_states][n_seqs]);

int main(int argc,char*argv[])
{
  //determine the number of MIP targets and the number of distinct paralogs
  FILE*miptargets=fopen(*(argv+1),"r");
	long num_mip_targets=0,num_plogs=0;
	count_targets(miptargets,&num_mip_targets,&num_plogs);

	//read in MIP specificities and target coordinates
	long coords[num_mip_targets];
	char specs[num_mip_targets][num_plogs+1];
	get_mip_info(miptargets,num_mip_targets,num_plogs+1,coords,specs);

	//set up vectors of copy number states and corresponding prior probabilities
	long para_cstates=strtol(*(argv+3),NULL,10)+1; //number of possible copy number states for each paralog
	long num_cstates=pow(para_cstates,num_plogs);
	double priors[num_cstates];
	long copy_states[num_cstates][num_plogs];
	set_priors(num_plogs,para_cstates,num_cstates,priors,copy_states);

	//allocate storage for likelihood graph
	struct node*likelihood_graph;
  likelihood_graph=(struct node*)malloc(num_mip_targets*num_cstates*sizeof(struct node));

	//set up output files
	char basename[101];
	FILE*outputs[3];
	strncpy(basename,*(argv+2),88);	
	init_output(outputs,basename,num_plogs);

	//for each individual, read in paralog-specific MIP read counts, calculate and store individual likelihoods of data for each MIP
	//under each possible copy number state, and use dynamic programming to infer paralog-specific copy number genotypes
	FILE*mipcounts=fopen(*(argv+2),"r");
	char indiv[101];
	double max0,nextmax0,max1,max2;
	long imax0,imax1,imax2,indiv_states[3],tmips[2];
	int num_trans;
	while(getc(mipcounts)!='\n')
		continue;
	while(fscanf(mipcounts,"%s",indiv)==1) //while there remains data to process
	{
		fseek(mipcounts,(-1*strlen(indiv)),SEEK_CUR);
		
		//initialize likelihood graph counts and log-likelihoods for a new individual
		init_graph(likelihood_graph,num_mip_targets,num_cstates);
		
		//process MIP data: get all relevant counts and calculate log-likelihoods of observed counts at each MIP under each possible copy number state
		fill_graph(likelihood_graph,mipcounts,num_plogs,num_mip_targets,num_cstates,num_plogs+1,specs,copy_states);

		//use dynamic programming to compute log-likelihoods for each copy number state and for best 1-transition and 2-transition paths ending at each copy number state
		parse_graph(likelihood_graph,num_plogs,num_mip_targets,num_cstates,priors,copy_states);

		//determine maximally likely paths having 0, 1, and 2 copy number state transitions and their corresponding log-likelihoods
		get_path_maxes(likelihood_graph,num_mip_targets,num_cstates,&max0,&nextmax0,&max1,&max2,&imax0,&imax1,&imax2);	

		//evaulate the evidence for multiple copy number states
		num_trans=imax(assess_path(likelihood_graph,max0,max1,max2,imax1,imax2,num_mip_targets,num_cstates,indiv_states,tmips,2),assess_path(likelihood_graph,max0,max1,max2,imax1,imax2,num_mip_targets,num_cstates,indiv_states,tmips,1));

		//print genotype information for the individual just analyzed
		print_output(outputs,indiv,max0,nextmax0,max1,max2,imax0,num_trans,indiv_states,tmips,num_mip_targets,coords,num_cstates,num_plogs,copy_states);
	}

	//clean up and exit
	free(likelihood_graph);
	fclose(miptargets);
	fclose(mipcounts);
	fclose(outputs[0]);
	fclose(outputs[1]);
	fclose(outputs[2]);	
	return 0;
}

//functions below
void count_targets(FILE*miptargs,long*nummips,long*numseqs)
{
	char spec[51];
	fpos_t pos;
	fgetpos(miptargs,&pos);
	while(fscanf(miptargs,"%*s %*s %*s %*s %*s %*s %s %*s %*s",spec)==1)
  {
  	(*nummips)++;
  }
	*numseqs=strlen(spec);
	fsetpos(miptargs,&pos);
	return;
}

void get_mip_info(FILE*miptargs,long nummips,long speclength,long targlocs[nummips],char specvecs[nummips][speclength])
{
	long i,start,end;
	for(i=0;i<nummips;i++)
	{
		fscanf(miptargs,"%*s %ld %*c %ld %*s %*s %*s %*s %s %*s %*s",&start,&end,specvecs[i]);
		targlocs[i]=(start+end)/2;
	}
	return;
}

void set_priors(long numseqs,long parastates,long numstates,double priorvec[numstates],long cstates[numstates][numseqs])
{
	long i,s,num;
	double pprobs[parastates];
	for(i=0;i<parastates;i++)
	{
		pprobs[i]=exp(L_PRIOR*labs(2-i)); //vector of prior probabilities for one paralog having each copy number state (assume CN = 2 is most likely)
	}
	for(i=0;i<numstates;i++)
	{
		num=i;
		priorvec[i]=0;
		for(s=(numseqs-1);s>=0;s--)
		{
			cstates[i][s]=num%parastates;
			priorvec[i]+=log(pprobs[num%parastates]);
			num/=parastates;
		}
	}
	return;
}

void init_output(FILE*outfiles[3],char*base,long numseqs)
{
	char*extensions[3]={".cncalls",".compevents",".simplecalls"};
	long i;
	for(i=0;i<3;i++)
	{
		base[strrchr(base,'.')-base]='\0';
		strcat(base,extensions[i]);
		outfiles[i]=fopen(base,"w");
	}
	fprintf(outfiles[2],"Individual\t");
	for(i=0;i<numseqs;i++)
	{
		fprintf(outfiles[2],"Para%ld_CN\t",i+1);
	}	
	fprintf(outfiles[2],"LOD_Score\tPossible_Complex_CN_Genotype\n");
	return;
}

void init_graph(struct node*graph,long nummips,long numstates)
{
	long i;
	for(i=0;i<(nummips*numstates);i++)
	{
		graph[i].likelihood=L_MIN;
		graph[i].max_0t=L_MIN;
		graph[i].max_1t=L_MIN;
		graph[i].max_2t=L_MIN;
		graph[i].path_1t=NULL;
		graph[i].path_2t=NULL;			
	}
	return;
}

void fill_graph(struct node*graph,FILE*mcounts,long numseqs,long nummips,long numstates,long speclength,char specvecs[nummips][speclength],long cnstates[numstates][numseqs])
{
	long mip,cstate,seq,copynum;
	unsigned int rawcounts[numseqs],counts[numseqs];
	double probs[numseqs];
	for(mip=0;mip<nummips;mip++)
	{
		//initialize counts
		init_counts(numseqs,rawcounts,counts);
		
		//read in raw counts and convert to final counts based on the specificity of the MIP
		fscanf(mcounts,"%*s %*s %*s %*s");
		for(seq=0;seq<numseqs;seq++)
		{
			fscanf(mcounts,"%u",&(rawcounts[seq]));
			counts[specvecs[mip][seq]-'A']+=rawcounts[seq]; //consolidate counts from identical sequences
		}
		
		//calculate likelihoods of observed counts at each MIP under each possible copy number state
		for(cstate=0;cstate<numstates;cstate++)
		{
			init_probs(numseqs,probs);
			copynum=num_copies(cstate,numseqs,numstates,cnstates);
			for(seq=0;seq<numseqs;seq++)
			{
				probs[specvecs[mip][seq]-'A']+=(double)cnstates[cstate][seq]/(double)copynum;
			}
			graph[mip*numstates+cstate].likelihood=gsl_ran_multinomial_lnpdf(numseqs,probs,counts);
			if(graph[mip*numstates+cstate].likelihood<L_MIN)
			{
				graph[mip*numstates+cstate].likelihood=L_MIN;
			}
		}
	}
	return;
}

void init_counts(long nseqs,unsigned int rcounts[nseqs],unsigned int fcounts[nseqs])
{
	long i;
	for(i=0;i<nseqs;i++)
	{
		rcounts[i]=0;
		fcounts[i]=0;
	}
	return;
}

long num_copies(long state,long nseqs,long nstates,long copystates[nstates][nseqs])
{
	long copynum=0,s;
	for(s=0;s<nseqs;s++)
	{
		copynum+=copystates[state][s];
	}
	return copynum;
}

void init_probs(long nseqs,double pvec[nseqs])
{
  long i;
  for(i=0;i<nseqs;i++)
  {
    pvec[i]=0.0;
  }
  return;
}

void parse_graph(struct node*graph,long numseqs,long nummips,long numstates,double priorvec[numstates],long cnstates[numstates][numseqs])
{
	long i;
	for(i=0;i<(nummips*numstates);i++)
	{
		if(!(i/numstates))
		{
			graph[i].max_0t=priorvec[i]+graph[i].likelihood; //apply priors to likelihoods for each copy number state at first MIP
			graph[i].max_1t=-DBL_MAX;
			graph[i].max_2t=-DBL_MAX;
		}
		else
		{
			graph[i].max_0t=graph[i].likelihood+graph[i-numstates].max_0t;
			graph[i].max_1t=graph[i].likelihood+segmax(graph,i,numseqs,numstates,cnstates,1);
			graph[i].max_2t=graph[i].likelihood+segmax(graph,i,numseqs,numstates,cnstates,2);
			if((i/numstates)==(nummips-1)) //for genotypes with at least one copy number state transition, apply priors to end as well
    	{
      	graph[i].max_1t+=priorvec[i%numstates];
      	graph[i].max_2t+=priorvec[i%numstates];
    	}		
		}
	}
	return;
}

double segmax(struct node*lgraph,long index,long nseqs,long nstates,long copystates[nstates][nseqs],int path)
{
	double max,tmax;
	struct node*prev;
	long j;
	max=((path==1)?lgraph[index-nstates].max_1t:lgraph[index-nstates].max_2t);
	prev=&(lgraph[index-nstates]);
	for(j=(index-nstates-index%nstates);j<(index-index%nstates);j++)
	{
		tmax=((path==1)?lgraph[j].max_0t:lgraph[j].max_1t);
		if((tmax>max)&&(trans_good(index,j,nstates,nseqs,copystates)))
		{
			max=tmax;
			prev=&(lgraph[j]);
		}
	}
	if(path==1) lgraph[index].path_1t=prev; else lgraph[index].path_2t=prev;
	return max;
}

int trans_good(long index,long index2,long n_states,long n_seqs,long pscn_states[n_states][n_seqs])
{
	long changes=0,copies=0,copies2=0,k;
	for(k=0;k<n_seqs;k++)
	{
		if(pscn_states[index%n_states][k]!=pscn_states[index2%n_states][k])
		{
			changes++;
		}
		copies+=pscn_states[index%n_states][k];
		copies2+=pscn_states[index2%n_states][k];
	}
	return ((changes==1)||((changes==2)&&(copies==copies2))); //transition is valid if paralog-specific copy number genotype change is due to deletion or duplication of a single paralog or interlocus gene conversion
}

void get_path_maxes(struct node*lgraph,long nmips,long nstates,double*m0,double*nm0,double*m1,double*m2,long*im0,long*im1,long*im2)
{
	long i;
	*m0=-DBL_MAX;
	*nm0=-DBL_MAX;
	*m1=-DBL_MAX;
	*m2=-DBL_MAX;
	*im0=(nmips*nstates-nstates);
	*im1=*im0;
	*im2=*im0;
	for(i=(nmips*nstates-nstates);i<(nmips*nstates);i++)
	{
		if(lgraph[i].max_0t>*m0)
		{
			if(diff_support(lgraph,i,*im0,nstates)) //ensure data support for the two most likely copy number states assessed so far differs (otherwise difference in their likelihoods is due to priors alone)
			{
				*nm0=*m0;
			}
			*m0=lgraph[i].max_0t;
			*im0=i;
		}
		else if(lgraph[i].max_0t>*nm0)
		{
			if(diff_support(lgraph,i,*im0,nstates))
			{
				*nm0=lgraph[i].max_0t;
			}
		}
		if(lgraph[i].max_1t>*m1)
		{
			*m1=lgraph[i].max_1t;
			*im1=i;
		}
		if(lgraph[i].max_2t>*m2)
    {
      *m2=lgraph[i].max_2t;
      *im2=i;
    }
	}
	return;
}

int diff_support(struct node*l_graph,long index,long index2,long n_states)
{
	int diff=0;
	long j=index,k=index2;
	for(j=index;j>=0;j-=n_states)
	{
		if(l_graph[j].likelihood!=l_graph[k].likelihood)
		{
			diff++;
			break;
		}
		k-=n_states;
	}
	return diff;
}

int assess_path(struct node*lgraph,double m0,double m1,double m2,long im1,long im2,long nmips,long nstates,long*cnstates,long*edgemips,int path)
{
	double ldiff=((path==1)?(m1-m0):(m2-dmax(m0,m1)));
	long statemips=0,numtrans=0,minmips=nmips,states[3]={-1,-1,-1},emips[2]={-1,-1},newstate;
	struct node*mip=((path==1)?&(lgraph[im1]):&(lgraph[im2]));	
	while(path)
	{
		statemips++;
		states[path]=(mip-&(lgraph[0]))%nstates;
		emips[path-1]=((mip-&(lgraph[0]))/nstates);
		mip=((path==1)?mip->path_1t:mip->path_2t);
		newstate=(mip-&(lgraph[0]))%nstates;
		if(states[path]!=newstate)
		{
			numtrans++;
			minmips=(statemips<minmips)?statemips:minmips;
			statemips=0;
			path--;
			if(!(path))
			{
				statemips=((mip-&(lgraph[0]))/nstates)+1;
				minmips=(statemips<minmips)?statemips:minmips;
				states[path]=(mip-&(lgraph[0]))%nstates;
			}
		}
	}
	if((ldiff>L_LRT)&&(minmips>=M_MIN)) //if there is evidence for multiple copy number states (see comment below for heuristics), store which states they are and which MIPs are located at transition points
	{
		edgemips[0]=emips[0]; //first MIP in second copy number state
		edgemips[1]=emips[1]-1; //last MIP in second copy number state
		cnstates[0]=states[0];
		cnstates[1]=states[1];
		cnstates[2]=states[2];
	}
	return (int)numtrans*(int)((ldiff>L_LRT)&&(minmips>=M_MIN)); //heuristics: A) likelihood difference between path assessed and each path with fewer transitions > 40, B) number of MIPs in each copy number state >= 5
}

int imax(int i1,int i2)
{
  return (i1>i2)?i1:i2;
}

double dmax(double d1,double d2)
{
  return (d1>d2)?d1:d2;
}

void print_output(FILE*outfiles[3],char*individual,double m0,double nm0,double m1,double m2,long im0,int ntrans,long*cnstates,long*edgemips,long nummips,long targlocs[nummips],long nstates,long nseqs,long copystates[nstates][nseqs])
{
	double lodscore=dmin((m0-nm0),LOD_MAX);
	long k;
	fprintf(outfiles[0],"Individual: %s\nL0=%lf L1-L0=%lf L2-L0=%lf\nGenotype (assuming no internal events): ",individual,m0,m1-m0,m2-m0);
	print_pscn(outfiles[0],im0%nstates,nstates,nseqs,copystates);
	fprintf(outfiles[0],"\n");
	if(ntrans)
	{
		printf(KYEL "\n%s may have a complex genotype!!! Graphically view the data for %s to ensure genotyping accuracy.\n",individual,individual); //alert user to possible complex genotype
		fprintf(outfiles[1],"%s\tL0=%lf L1-L0=%lf L2-L0=%lf\tGenotype (assuming no internal events): ",individual,m0,m1-m0,m2-m0); //print to 2nd output file so you can easily visualize only these individuals
		print_pscn(outfiles[1],im0%nstates,nstates,nseqs,copystates);
		fprintf(outfiles[1],"\t");
		fprintf(outfiles[2],"%s\t",individual);
		for(k=0;k<nseqs;k++)
    {
      fprintf(outfiles[2],"%ld\t",copystates[im0%nstates][k]);
    }	
    fprintf(outfiles[2],"%lf\tYES\n",lodscore);
		if(ntrans==1)
		{
			fprintf(outfiles[0],"Genotype (assuming 1 transition): ");
			print_pscn(outfiles[0],cnstates[0],nstates,nseqs,copystates);
			fprintf(outfiles[0]," ");
			print_pscn(outfiles[0],cnstates[1],nstates,nseqs,copystates);
			fprintf(outfiles[0],", First MIP in State 2: %ld\n\n",targlocs[edgemips[0]]);
			fprintf(outfiles[1],"Genotype (assuming 1 transition): ");
      print_pscn(outfiles[1],cnstates[0],nstates,nseqs,copystates);
      fprintf(outfiles[1]," ");
      print_pscn(outfiles[1],cnstates[1],nstates,nseqs,copystates);
      fprintf(outfiles[1],", First MIP in State 2: %ld\n",targlocs[edgemips[0]]);
		}
		else
		{
			fprintf(outfiles[0],"Genotype (assuming 2 transitions): ");
      print_pscn(outfiles[0],cnstates[0],nstates,nseqs,copystates);
      fprintf(outfiles[0]," ");
      print_pscn(outfiles[0],cnstates[1],nstates,nseqs,copystates);
			fprintf(outfiles[0]," ");
			print_pscn(outfiles[0],cnstates[2],nstates,nseqs,copystates);
      fprintf(outfiles[0],", First MIP in State 2: %ld, Last MIP in State 2: %ld\n\n",targlocs[edgemips[0]],targlocs[edgemips[1]]);
			fprintf(outfiles[1],"Genotype (assuming 2 transitions): ");
      print_pscn(outfiles[1],cnstates[0],nstates,nseqs,copystates);
      fprintf(outfiles[1]," ");
      print_pscn(outfiles[1],cnstates[1],nstates,nseqs,copystates);
      fprintf(outfiles[1]," ");
      print_pscn(outfiles[1],cnstates[2],nstates,nseqs,copystates);
      fprintf(outfiles[1],", First MIP in State 2: %ld, Last MIP in State 2: %ld\n",targlocs[edgemips[0]],targlocs[edgemips[1]]);
		}
	}
	else
	{
		fprintf(outfiles[0],"\n");
		fprintf(outfiles[2],"%s\t",individual);
		for(k=0;k<nseqs;k++)
		{
			fprintf(outfiles[2],"%ld\t",copystates[im0%nstates][k]);
		}
		fprintf(outfiles[2],"%lf\tNO\n",lodscore);
	}
	return;
}

double dmin(double d1,double d2)
{
  return (d1<d2)?d1:d2;
}

void print_pscn(FILE*out,long state,long n_states,long n_seqs,long pscn_states[n_states][n_seqs])
{
	long s;
	for(s=0;s<n_seqs;s++)
  {
    fprintf(out,"%ld",pscn_states[state][s]);
  }
	return;	
}
