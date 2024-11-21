/******************************************************************************/
/* PROD PROGRAM:    derived_rna.sas
/* WORK PROGRAM:    /home/phacs/HOPE/programs/Yoshi/Programs/derived_rna.sas
/* 
/* RISK ASSESSMENT: LOW
/* DOCUMENT ID:     PROG.10083.v4
/* 
/* PURPOSE:         To create a comprehensive RNA profile for each HOPE participant 
                    by supplementing RNA data from HOPE with RNA data from other 
                    protocols for HOPE participants co-enrolled in other studies. 
/* 
/* SOURCE PRGM:     none
/* INPUT:           derv100.matrna, derv300.comb_rnamrg, perm400.hxw0219, permhope.hxw0233
                    permhope.anstab, permhope.QLW0360, permhope.evw0379, permhope.evw0380
                    derv100.master, dervhope.master
/* OUTPUT:          out.rna
/* MACROS USED:     none
/* EXEMPTIONS:      none
/* 
/* AUTHOR:          Yashodhan (Yoshi) Aher
/* CREATION DATE:   05/28/24
/* 
/* NOTES:           none
/* MODIFICATIONS:   none
/******************************************************************************/

/*****************************************************************************************************************************************************************************
                                                                        REVIEW OF RELEVANT INPUT DATASETS
*****************************************************************************************************************************************************************************/

title 'SMARTT Master DD'; 
proc contents data=derv100.master;
run; 

title 'SMARTT RNA DD'; 
proc contents data=derv100.matrna; 
run; 

title 'AMP Up RNA DD'; 
proc contents data=derv300.comb_rnamrg; 
run; 

title 'AMP Up Lite Viral Load CRF'; 
proc contents data=perm400.hxw0219; 
run;

title 'HOPE Master DD'; 
proc contents data=dervhope.master; 
run; 

title 'HOPE RNA CRF Data';
proc contents data=permhope.hxw0233;
run; 

title 'HOPE Form EVW0379';
proc contents data=permhope.evw0379; 
run; 

title 'HOPE Form EVW0380';
proc contents data=permhope.evw0380; 
run;

title 'HOPE Pregnancy History Questionnaire';
proc contents data=permhope.QLW0360; 
run;

title 'HOPE QLW0357';
proc contents data=permhope.QLW0357; 
run;

title 'HOPE Anstab';
proc contents data=permhope.anstab; 
run;

title 'HOPE Visit Form';
proc contents data=permhope.ph5805;  
run;


/*******************************************************************************************************************************
                                            AGGREGATING AND STANDARDIZING RNA DATA FROM ACROSS PROTOCOLS
********************************************************************************************************************************/

/*                                                    DIFFERENT QUANTCD FORMATTING
Study   Var         Type    Length  Format
HOPE    quantcd     Num     3       QUANT9FT.
SMARTT  qntcd       Num     8       QNTPH200F.
AMP Up  qntcd       Num     8       QNTPH200F.
AUL     quantfcd    Num     3       QUANT2FT.

value QNTPH200F                                                                                                                                                      
                  1='<' 
                  1.5 = 'Mean of 2 values with different quantifiers (</=)'                                                   
                  2='='                                                                                                   
                  3='>' 

 CBAR     QUANT9FT     
                        1      Equal to (=)
                        2      Greater than (>)
                        3      Less than (<)
                        4      Unknown


CBAR QUANT2FT.
                        1      Equal to (=)
                        2      Greater than (>)
                        3      Less than (<)
                        4      Unknown


*/


/*Quantcd has various formats across protocols so this is a standard format than ensures no data is lost across formats*/
proc format;
    value stan_qnt
        1 = 'Equal to (=)'
        2 = 'Greater than (>)'
        3 = 'Less than (<)'
        4 = 'Unknown'
        1.5 = 'Mean of 2 values with different quantifiers (</=)'
    ;

run; 

/*Keeping unique ids from HOPE master dataset because we only care about RNA data relevant to HOPE patids*/
proc sort data=dervhope.master out=master(keep=patid) nodupkey; 
    by patid; 
run; 

/*Get patid and biomomid from SMARTT master. This is needed because patid is recorded under childid in SMARTT, so 
  we cannot successfully bring in SMARTT data unless we merge HOPE and SMARTT using biomomid*/
proc sort data=derv100.master out=biomom(keep=patid biomomid where=(biomomid ne .));
    by patid biomomid;  
run; 

/*Creating variables found in HOPE, AUL, and SMARTT that were not already in AMP Up(comb_rnamrg)*/
data ampup_rna(drop=qntcd);  
  set derv300.comb_rnamrg(keep=patid copyml qntcd source specdt); 
  length assaysp $70 assaytyp 3 quantcd 3; 
  assaysp=strip(''); 
  assaytyp=.;
  source=strip('PH300'); 
  /*qntcd from ampup is formatted using QNTPH200F and HOPE uses QUANT9FT. so new var 'quantcd' reassigns 
    values based on QUANT9FT to ensure the data is consistent across studies in the final dataset*/
  if qntcd=1 then quantcd=3; 
    else if qntcd=2 then quantcd=1;
    else if qntcd=3 then quantcd=2;
    else if qntcd=1.5 then quantcd=1.5; 
run; 

/*Removing duplicate observations using following by-groups*/
proc sort data=ampup_rna nodupkey;  
    by patid specdt quantcd copyml; 
run;

/*Creating variables found in SMARTT, AMP Up, and HOPE that were not already in AUL viral load CRF*/   
data aul_rna; 
    set perm400.hxw0219(keep=assaydt assaysp assaytyp copyml patid quantfcd rename=(assaydt=specdt quantfcd=quantcd)); 
    length source $200 quantcd 3; 
    source=strip('PH400'); 
    /*quantcd doesn't need to be reassigned values because the format AUL uses, although different in name, is 
      actually the same as HOPE*/
run;   

/*Removing duplicate observations using following by-groups*/
proc sort data=aul_rna nodupkey;  
    by patid specdt quantcd copyml; 
run; 

/*Creating variables found in HOPE, AMP Up, and AUL that were not already in SMARTT RNA derived dataset*/
data smarttrna(drop=qntcd);
  set derv100.matrna(keep=patid copyml qntcd source specdt);
  length assaysp $70 assaytyp quantcd 3;
  assaysp=strip(''); 
  assaytyp=.;
  source=strip('PH100'); 
  /*qntcd from SMARTT is formatted using QNTPH200F and HOPE uses QUANT9FT. so new var 'quantcd' reassigns 
  values based on QUANT9FT to ensure the data is consistent across studies in the final dataset*/
  if qntcd=1 then quantcd=3; 
    else if qntcd=2 then quantcd=1;
    else if qntcd=3 then quantcd=2;
    else if qntcd=1.5 then quantcd=1.5; 
run;

/*Add biomomid to smartt_rna and rename to patid in order to merge dataset later
  (Reminder: SMARTT is recorded under child PIDs & HOPE under the mother's PID so
  if SMARTT PID is used there will be lots of missing data)*/
data smartt_rna(drop=patid rename=(biomomid=patid)); 
    merge smarttrna biomom; 
    by patid; 
run; 

/*Removing duplicate observations using following by-groups*/
proc sort data=smartt_rna nodupkey; 
    by patid specdt quantcd copyml; 
run; 

/*Creating variables found in SMARTT, AMP Up, and AUL that were not already in HOPE viral load CRF */
data hope_rna; 
    set permhope.hxw0233(keep=patid copyml quantcd specdt assaysp assaytyp); 
    length source $200; 
    source='HOPE'; 
run; 

/*Removing duplicate observations using following by groups*/
proc sort data=hope_rna nodupkey; 
    by patid specdt quantcd copyml; 
run; 

/*******************************************************************************************************************************
                                                            DATES PROGRAMMING
*******************************************************************************************************************************/

/*The mkvars macro essentially works like a transpose where if you had a dataset with 2 variables
  patid and pregdt and the dataset had 4 observations, 2 obs where patid=1 and 2 obs where patid = 2 then all the dates related
  to a patid are compuiled onto 1 line and if the last.patid statement is used it keeps the observation that has all the 
  pregdts relevant to that patid which makes comparison of dates possible. You will see this macro many times and it works
  the same way every time

  Original Data   What the macro does to the orig data  Final data when last.var is used

  patid pregdt   |   patid   pregdt_1    pregdt_2     |   patid   pregdt_1    pregdt_2
    1     2011   |     1         2011         .       |     1         2011        2012
    1     2012   |     1         2011         2012    |     2         2013        2014
    2     2013   |     2         2013         .       |        
    2     2014   |     2         2013         2014    |            
  
*/
%macro mkvars(counter);
    if counter=&counter then do; 
        anstab_pregdt_&counter=anstab_pregdt; anstab_gesage_&counter=anstab_gesage;   
    end; 
%mend mkvars;

%macro anstab(); 

/*ANSTAB_pregdt: Select the observations that have gestational ages */
data anstab_gesage;
    set permhope.anstab(where=(qnum=6 and resp1 ne '')); /*qnum=6 corresponds to the following question:
                                                           'Enter the participant's gestational age (in weeks) by best obstetric estimates.|*/
    anstab_gesage=input(resp1, best3.); 
    label anstab_gesage='Gestational Age from Anstab';
    keep patid anstab_gesage; 
run; 

/*ANSTAB_pregdt: Select the observations that have preg dates */
data anstab_pregdt;
    set permhope.anstab(where=(qnum=11 and resp1 ne '')); /*qnum=11 corresponds to the following question:
                                                           'Enter the participant's most recent delivery date.*/
    anstab_pregdt=input(resp1, anydtdte10.); 
    format anstab_pregdt date9.;
    label anstab_pregdt='Pregnancy Outcome Date from Anstab';
    keep patid anstab_pregdt; 
run; 

/*ANSTAB: Contains patid and either a gesage or pregdt for every patid in Anstab*/
data anstab; 
    merge anstab_gesage(in=ingesage) anstab_pregdt(in=inpregdt) dervhope.master(in=inmaster keep=patid stratcat2 randdt); 
    by patid;  
    if ingesage or inpregdt; 

    /*Subset for pregnant participants and for recently delivered or within 1-year postpartum participants. For the recently delivered or within 
      1-year postpartum participants, we want the "most recent" pregnancy which refers to the pregnancy with an outcome date closest to but not 
      after the randdt. 

      This block of code doesn't need to use entryvisdt/enrtystratcat2(which you will see used in QLW0360 and EVW0179) because anstab is tied to 
      enrollment(compared to the other sources that are based on visit). The example below will help to visualize this a bit better: 

      A participant enrolls as pregnant. This means in anstab, the participant will have gestational week but no pregdt (as the outcome didnt happen 
      as of randdt). The participants entry visit doesnt happen until the participant gives birth and is 3 months postpartum. Therefore, the 
      entrystratcat will be postpartum at the entryvisdt, and this participant will have data in the other forms that indicate pregdt, gestational age, 
      etc. to calculate the lmpdt. So effectively the info from anstab wont get used for this participant.

    */
    
    if stratcat2 = 2 or (stratcat2 in (3,4) and anstab_pregdt ne . and anstab_pregdt<=randdt);


    keep patid anstab_pregdt anstab_gesage; 

run; 

/*Store max number of pregnancies in a macro variable to be used later*/
proc sql noprint;
    select max(count) into :max_anstab_dts
    from (
        select patid, count(*) as count
        from anstab
        group by patid);
quit;


/*If a patid has more than 1 pregdt then further processing is needed to figure out which date to select*/
%if &max_anstab_dts > 1 %then %do; 
    
    proc sort data=anstab; 
        by patid; 
    run; 

    data anstab; 
        merge anstab(in=inanstab) dervhope.master(in=inmaster keep=patid stratcat2 randdt); 
        by patid;
        if inanstab;    

        /* Initialize new pregdt and gesage dummy variables */
        %do i = 1 %to &max_anstab_dts;
            length anstab_pregdt_&i 5 anstab_gesage_&i 3;
        %end; 

        /*Retain dummy pregdt and gesage vars to hold values of pregdts and gesage for a given patid*/
        retain anstab_pregdt_1-anstab_pregdt_%sysfunc(compress(&max_anstab_dts)) 
               anstab_gesage_1-anstab_gesage_%sysfunc(compress(&max_anstab_dts)) 
               counter; 
        
        /*Define arrays to store multiple pregdts and gesages*/ 
        array retarr_pregdts{*} anstab_pregdt_1-anstab_pregdt_%sysfunc(compress(&max_anstab_dts)); 
        array retarr_gesage{*} anstab_gesage_1-anstab_gesage_%sysfunc(compress(&max_anstab_dts));

        /*Call mkvars macro */
        %do i=1 %to &max_anstab_dts;
            %mkvars(&i); 
        %end; 

        drop i anstab_pregdt anstab_gesage; 

        /*Observation where last.patid is true is where all dummy pregdts and gesages will be found*/
        if last.patid;  
    run; 

    data anstab; 
        set anstab;
        /* Initialize the variables to assign final preg date and gesage to*/
        length anstab_pregdt 5 anstab_gesage 3; 
        anstab_pregdt = .;
        anstab_gesage = .;

        /* Define an array with the date variables */
        array pregdts{*} anstab_pregdt_1-anstab_pregdt_%sysfunc(compress(&max_anstab_dts));
        array gesages{*} anstab_gesage_1-anstab_gesage_%sysfunc(compress(&max_anstab_dts));
        
        /*
           This code identifies the pregnancy outcome date that is closest to but not after 
           the enrollment date (randdt) for participants who are recently delivered or within 
           1 year postpartum.
           
           - It initializes a variable `min_diff` to store the minimum difference between 
             pregnancy outcome dates and the enrollment date.
           
           - The code loops through an array of pregnancy dates (`pregdts`), checking each 
             non-missing date that is on or before the enrollment date (`randdt`).
           
           - For each valid date, the code calculates the absolute difference between the 
             pregnancy date and the enrollment date.
           
           - If the current difference is smaller than the previously recorded `min_diff` 
             (or if `min_diff` is missing, which happens on the first valid check), the 
             code updates `min_diff` and assigns the current pregnancy date to `anstab_pregdt`.
           
           - This ensures that `anstab_pregdt` holds the pregnancy date that is closest to 
             but not after the enrollment date.
        */

        min_diff = .;

        do i = 1 to dim(pregdts);
            if not missing(pregdts[i]) and pregdts[i] <= randdt then do; 
                diff = abs(randdt - pregdts[i]); 
                if missing(min_diff) or diff < min_diff then do; 
                    min_diff = diff;
                    anstab_pregdt = pregdts[i]; 
                end; 
            end;
        end;
 
        /* This code assigns the gestational age (`anstab_gesage`) that corresponds to the 
           pregnancy date (`anstab_pregdt`) found in the previous step.
           
           - The function `whichn` identifies the position (`x`) of `anstab_pregdt` in the 
             array of pregnancy dates (`pregdts[*]`).
           
           - If `anstab_pregdt` is not missing, the code uses the index `x` to assign the 
             corresponding gestational age from the `gesgaes` array to `anstab_gesage`. */

        x=whichn(anstab_pregdt, of pregdts[*]); 
        if anstab_pregdt ne . then anstab_gesage = gesgaes[x]; 

        keep patid anstab_pregdt anstab_gesage; 
    run; 

%end;
 
%mend anstab;  

%anstab(); 

%macro mkvars(counter);
    if counter=&counter then do; 
        pregdt0360_&counter=pregdt0360; gesage0360_&counter=gesage0360;   
    end; 
%mend mkvars;

%macro QLW0360(); 
/*Getting pregdt and gesage from QLW0360. */
data QLW0360; 
  merge dervhope.master(in=inmaster keep=patid entrystratcat2 stratcat2 randdt entryvisdt) 
        permhope.QLW0360(in=inqlw keep=patid pgoutdt1 pgdur1 pgoutdt2 pgdur2
                         where=(pgoutdt1 ne . or pgoutdt2 ne .)); 
  by patid; 
  if inqlw; 

  length pregdt0360 5 gesage0360 3;  
  
  pregdt0360 = .; /* Initialize pregdt0360 to hold the selected date*/

    /*
       This code selects the closest pregnancy outcome date (`pregdt0360`) on or before 
       the entry visit date (`entryvisdt`) from two possible dates: `pgoutdt1` and `pgoutdt2`.

       - The block runs only if `entrystratcat2` is not missing.
       
       - First, it checks if both `pgoutdt1` and `pgoutdt2` are non-missing.
         - If both dates are on or before the entry visit date, the code compares the 
           absolute differences between each date and `entryvisdt`. The date closer to 
           `entryvisdt` is assigned to `pregdt0360`.
         - If only one of the dates is on or before `entryvisdt`, that date is assigned 
           to `pregdt0360`.
       
       - If only one of `pgoutdt1` or `pgoutdt2` is non-missing and on or before 
         `entryvisdt`, that date is directly assigned to `pregdt0360`.
       
       - This logic ensures that `pregdt0360` is the pregnancy outcome date closest to 
         the entry visit date, provided it does not occur after the entry visit.
    */

  
  %if entrystratcat2 ne . %then %do; 
    if not missing(pgoutdt1) and not missing(pgoutdt2) then do;

        if pgoutdt1 <= entryvisdt and pgoutdt2 <= entryvisdt then do;
                
            if abs(entryvisdt - pgoutdt1) < abs(entryvisdt - pgoutdt2) then pregdt0360 = pgoutdt1;
            else pregdt0360 = pgoutdt2;

        end;
     
        else if pgoutdt1 <= entryvisdt then pregdt0360 = pgoutdt1;  
        else if pgoutdt2 <= entryvisdt then pregdt0360 = pgoutdt2;

    end;
    
    else if not missing(pgoutdt1) and pgoutdt1 <= entryvisdt then pregdt0360 = pgoutdt1;
    else if not missing(pgoutdt2) and pgoutdt2 <= entryvisdt then pregdt0360 = pgoutdt2;
  %end; 
  /*In the event that the entryvisit did not occur, the code below does the same as the block above except it compares the pgoutdt vars to the randdt*/
  %else %if entrystratcat2 = . %then %do;
    if not missing(pgoutdt1) and not missing(pgoutdt2) then do;

        if pgoutdt1 <= randdt and pgoutdt2 <= randdt then do;
                
            if abs(randdt - pgoutdt1) < abs(randdt - pgoutdt2) then pregdt0360 = pgoutdt1;
            else pregdt0360 = pgoutdt2;

        end;

        else if pgoutdt1 <= randdt then pregdt0360 = pgoutdt1;
   
        else if pgoutdt2 <= randdt then pregdt0360 = pgoutdt2;

    end;

    else if not missing(pgoutdt1) and pgoutdt1 <= randdt then pregdt0360 = pgoutdt1;
    else if not missing(pgoutdt2) and pgoutdt2 <= randdt then pregdt0360 = pgoutdt2;  
  %end; 
  
  /*Match the gesage with the corresponding pgoutdt*/
  if      pregdt0360=pgoutdt1 then pgdur=pgdur1; 
  else if pregdt0360=pgoutdt2 then pgdur=pgdur2; 

  /* Set gesage to larger approximate gesage(refer to CRF if needed) according to Jessica's 
     "Calculating threshold window for earliest RNA CD4" word doc found in: home/phacs/HOPE/documents/programming*/
  if      pgdur=1 then gesage0360=8; 
  else if pgdur=2 then gesage0360=13;
  else if pgdur=3 then gesage0360=26;
  else if pgdur=4 then gesage0360=36; 
  else if pgdur=5 then gesage0360=40;
  else if pgdur=6 then gesage0360=.; /*6 corresponds to 'Unknown' in the CRF*/

  if gesage0360 = . then gesage0360 = 40; /*Jessica's document says if gestational age is missing, assume 40 weeks gestation*/

  keep patid pregdt0360 gesage0360 entrystratcat2 stratcat2 randdt entryvisdt; 

  format pregdt0360 date9.;
  label pregdt0360='Pregnancy Outcome Date from QLW0360'
        gesage0360='Gestational Age from QLW0360'; 
run; 


/*Store max number of pregnancies in a macro variable to be used later*/
proc sql noprint;
    select max(count) into :max_qlw0360_dts
    from (
        select patid, count(*) as count
        from QLW0360
        group by patid
        having count(*) > 1
        );
quit;

/*If a participant had more than 2 pregancies within the last 12 months, a second sequence of this CRF was filled out which means
  there's atleast a third date of pregnancy outcome that has to be looked at so further processing is required*/
%if &max_qlw0360_dts > 1 %then %do; 
    
    proc sort data=qlw0360; 
        by patid; 
    run;

    data qlw0360; 
        set qlw0360; 
        by patid; 
        /* Create the new variables */
        %do i = 1 %to &max_qlw0360_dts;
            length pregdt0360_&i 5 gesage0360_&i 3;
        %end; 

        /*Retain dummy pregdt and gesage vars to hold values of pregdts and gesage for a given patid*/
        retain pregdt0360_1-pregdt0360_%sysfunc(compress(&max_qlw0360_dts)) 
               gesage0360_1-gesage0360_%sysfunc(compress(&max_qlw0360_dts)) 
               counter;    
        
        /*Define arrays to store multiple pregdts and gesages*/ 
        array pregdts{*} pregdt0360_1-pregdt0360_%sysfunc(compress(&max_qlw0360_dts)); 
        array gesages{*} gesage0360_1-gesage0360_%sysfunc(compress(&max_qlw0360_dts));

        /*Initialize dummy vars*/
        if first.patid then do; 
            counter=1;
            do i=1 to dim(pregdts); 
               pregdts(i)=.; 
               gesages(i)=.;
            end; 
        end; 
        else counter+1; 

        /*Call mkvars macro */
        %do i=1 %to &max_qlw0360_dts; 
            %mkvars(&i); 
        %end; 

        drop i pregdt0360 gesage0360; 

        /*Observation where last.patid is true is where all dummy pregdts and gesages will be found*/
        if last.patid;  
    run; 
    
    data qlw0360; 
        set qlw0360;
        /* Initialize the variable to assign final date and gesage to*/
        length pregdt0360 5 gesage0360 3; 
        pregdt0360 = .;
        gesage0360 = .;

        /* Define an array with the date variables */
        array pregdts[*] pregdt0360_1-pregdt0360_%sysfunc(compress(&max_qlw0360_dts));
        array gesages{*} gesage0360_1-gesage0360_%sysfunc(compress(&max_qlw0360_dts));
        
        /* Initialize the minimum difference */
        min_diff = .;

        /*
           This code identifies the closest pregnancy date (`pregdt0360`) on or before 
           the entry visit date (`entryvisdt`) from an array of pregnancy dates (`pregdts`).

           - The block runs only if `entrystratcat2` is not missing.
           
           - The code loops through each element in the array `pregdts`.
             - For each date, it checks if the date is non-missing and occurs on or before 
               the entry visit date (`entryvisdt`).
             - It calculates the absolute difference (`diff`) between the current date 
               and `entryvisdt`.
             - If this difference is smaller than the current minimum difference (`min_diff`) 
               (or if `min_diff` is missing, indicating the first valid date), the code 
               updates `min_diff` and assigns the current date to `pregdt0360`.

           - This process ensures that `pregdt0360` holds the date from the array that is 
             closest to but not after the entry visit date.
        */


        %if entrystratcat2 ne . %then %do; 
            do i = 1 to dim(pregdts);
                if not missing(pregdts[i]) and pregdts[i] <= entryvisdt then do; 
                    diff = abs(entryvisdt - pregdts[i]); 
                    if missing(min_diff) or diff < min_diff then do; 
                        min_diff = diff;
                        pregdt0360 = pregdts[i]; 
                    end; 
                end;
            end;
        %end; 
        /*If entryvisit did not happen, then the following block of code does the same as above except pregdts are compared to randdt instead of entryvisitdt*/
        %else %if entrystratcat2 = . %then %do; 
            do i = 1 to dim(pregdts);
                if not missing(pregdts[i]) and pregdts[i] <= randdt then do; 
                    diff = abs(randdt - pregdts[i]); 
                    if missing(min_diff) or diff < min_diff then do; 
                        /*If so, it updates min_diff and sets pregdt0360 to the current date*/
                        min_diff = diff;
                        pregdt0360 = pregdts[i]; 
                    end; 
                end;
            end;
        %end; 
        
        /* This code assigns the gestational age (`gesage0360`) corresponding to the 
           closest pregnancy date (`pregdt0360`) from the array of pregnancy dates (`pregdts`).

           - The `whichn` function finds the position (`x`) of `pregdt0360` within the array 
             `pregdts[*]`.
           
           - If `pregdt0360` is not missing, the code uses the index `x` to assign the 
             corresponding gestational age from the `gesages` array to `gesage0360`. */

        x=whichn(pregdt0360, of pregdts[*]); 

        if pregdt0360 ne . then gesage0360 = gesages[x]; 

        /*Jessica's document says if gestational age is missing, assume 40 weeks gestation*/
        if gesage0360 = . then gesage0360 = 40; 

        keep patid pregdt0360 gesage0360; 
    run;     

%end;  

%mend QLW0360; 

%qlw0360();

%macro mkvars(counter);
    if counter=&counter then do; 
        pregdt0379_&counter=pregdt0379; gesage0379_&counter=gesage0379;   
    end; 
%mend mkvars;

%macro evw0379(); 

proc sql noprint;
    /*Store all pregdt1, pregdt2, etc. vars in a macro variable called pregdt_vars*/    
    select name 
    into :pregdt_vars separated by ' '
    from dictionary.columns
    where libname='PERMHOPE' 
      and memname=upcase("EVW0379") 
      and upcase(name) like upcase("pregdt%")
      and upcase(name) ne upcase("pregdt");

    /* Count the number of pregdt variables */
    select count(name) 
    into :num_pregdt_vars
    from dictionary.columns
    where libname='PERMHOPE' 
      and memname=upcase("EVW0379") 
      and upcase(name) like upcase("pregdt%")
      and upcase(name) ne upcase("pregdt");  

    /*Store all f1gesage, f2gesage, etc. vars in a macro variable called gesage_vars*/  
    select name 
    into :gesage_vars separated by ' '
    from dictionary.columns
    where libname='PERMHOPE' 
      and memname=upcase("EVW0379") 
      and upcase(name) like '%GESAGE%';

quit; 

/*Getting pregdt and gesage from evw0379. */
data evw0379;
  merge dervhope.master  (in=inmaster keep=patid entrystratcat2 stratcat2 randdt entryvisdt) 
        permhope.evw0379 (in=inevw keep=patid entryvis &pregdt_vars &gesage_vars); 
   by patid;

   if inevw;

   /* Filter for obs where atleast one pregdt is not missing */
   if cmiss(of &pregdt_vars.) < &num_pregdt_vars.;

   /*This CRF is also completed at follow-up visits but for date calculation purposes you only need the entry visit form(s)*/
   if entryvis = 1; 

   length pregdt0379 5 gesage0379 3; 

   /*Define arrays to store date of pregnancy outcome for each fetus and corresponding gesages*/
   array pregdts {*} &pregdt_vars;
   array gesages {*} &gesage_vars;

   /*The `coalesce` function selects the first non-missing value from the list 
     of pregnancy dates (`pregdt1` to `pregdt<num_pregdt_vars>`), dynamically 
     determined by the macro variable `&num_pregdt_vars`; these are potentially 
     twins/triplets, etc. so since the dates will be nearly identical, it
     doesn't matter which one we use.
       
     The `whichn` function finds the position (`x`) of `pregdt0379` within the 
     `pregdts` array.
       
     If `pregdt0379` is not missing, the code assigns the corresponding 
     gestational age from the `gesages` array to `gesage0379` using the index `x`.*/

   pregdt0379 = coalesce(of pregdt1-pregdt%sysfunc(compress(&num_pregdt_vars)));
   x=whichn(pregdt0379, of pregdts[*]); 
   if pregdt0379 ne . then gesage0379=gesages[x]; 

   /*Subset the data based on if pregdt0379 is <= entryvisitdt or <= randdt if entryvisitdt is missing*/
   %if entrystratcat2 ne . %then %do; 
      if pregdt0379<=entryvisdt; 
   %end; 
   %else %if entrystratcat2 = . %then %do; 
      if pregdt0379<=randdt;
   %end; 

   keep patid stratcat2 randdt entryvisdt pregdt0379 gesage0379;

   format pregdt0379 date9.;
run; 

/*Store max number of pregnancies in a macro variable to be used later*/
proc sql;
    select max(count) into :max_evw0379_dts
    from (
        select patid, count(*) as count
        from evw0379
        group by patid
        );
quit;

/*If more than 1 pregnancy for a patid, then there's more processing required to evaluate all pregnancy outcome dates for a given patid*/
%if &max_evw0379_dts > 1 %then %do; 
    
    proc sort data=evw0379; 
        by patid; 
    run; 

    data evw0379; 
        set evw0379; 
        by patid; 
        /* Create the new variables */
        %do i = 1 %to &max_evw0379_dts;
            length pregdt0379_&i 5 gesage0379_&i 3;
        %end; 

        /*Retain dummy pregdt and gesage vars to hold values of pregdts and gesage for a given patid*/
        retain pregdt0379_1-pregdt0379_%sysfunc(compress(&max_evw0379_dts)) 
               gesage0379_1-gesage0379_%sysfunc(compress(&max_evw0379_dts)) 
               counter;    
        
        /*Define arrays to store multiple pregdts and gesages*/ 
        array pregdts{*} pregdt0379_1-pregdt0379_%sysfunc(compress(&max_evw0379_dts)); 
        array gesages{*} gesage0379_1-gesage0379_%sysfunc(compress(&max_evw0379_dts));

        if first.patid then do; 
            counter=1;
            do i=1 to dim(pregdts); 
               pregdts(i)=.; 
               gesages(i)=.;
            end; 
        end; 
        else counter+1;  

        /*Call mkvars macro */
        %do i=1 %to &max_evw0379_dts; 
            %mkvars(&i); 
        %end; 

        drop i pregdt0379 gesage0379; 

        /*Observation where last.patid is true is where all dummy pregdts and gesages will be found*/
        if last.patid; 
    run; 

    data evw0379; 
        set evw0379;
        /* Initialize the variable to assign final date and gesage to*/
        length pregdt0379 5 gesage0379 3; 
        pregdt0379 = .;
        gesage0379 = .;
        /* Define an array with the date variables */
        array pregdts[*] pregdt0379_1-pregdt0379_%sysfunc(compress(&max_evw0379_dts));
        array gesages{*} gesage0379_1-gesage0379_%sysfunc(compress(&max_evw0379_dts));
        
        /* Initialize the minimum difference */
        min_diff = .;
        
        /*This code block checks if `entrystratcat2` is not missing and, if so, 
           iterates through the array of pregnancy dates (`pregdts`) to find the 
           closest date that is on or before `entryvisdt`.

           - For each date, it verifies that the date is non-missing and less than 
             or equal to `entryvisdt`, calculating the absolute difference between 
             `entryvisdt` and the current date.

           - It checks if this difference is less than the current minimum difference 
             (`min_diff`) or if `min_diff` is missing, indicating the first valid date.

           - If the current date is closer, it updates `min_diff` and assigns the 
             current pregnancy date to `pregdt0379`.*/


        %if entrystratcat2 ne . %then %do; 
            do i = 1 to dim(pregdts);
                if not missing(pregdts[i]) and pregdts[i] <= entryvisdt then do; 
                    diff = abs(entryvisdt - pregdts[i]); 
                    if missing(min_diff) or diff < min_diff then do; 
                        /*If so, it updates min_diff and sets pregdt0379 to the current date*/
                        min_diff = diff;
                        pregdt0379 = pregdts[i]; 
                    end; 
                end;
            end;
        %end;
        /*If entryvisit did not happen, then the following block of code does the same thing as the above except it compares pregdts to randdt instead of entryvisitdt*/
        %else %if entrystratcat2 = . %then %do; 
            do i = 1 to dim(pregdts);
                if not missing(pregdts[i]) and pregdts[i] <= randdt then do; 
                    diff = abs(randdt - pregdts[i]); 
                    if missing(min_diff) or diff < min_diff then do; 
                        min_diff = diff;
                        pregdt0379 = pregdts[i]; 
                    end; 
                end;
            end;
        %end; 
        
        /* This code segment determines the index of the closest valid pregnancy date 
           (`pregdt0379`) within the array of pregnancy dates (`pregdts`).

           - The `whichn` function is used to find the position (`x`) of `pregdt0379` 
             in the `pregdts` array.

           - If `pregdt0379` is not missing, the code assigns the corresponding 
             gestational age from the `gesages` array to `gesage0379` using the index `x`.*/

        x=whichn(pregdt0379, of pregdts[*]); 
        if pregdt0379 ne . then gesage0379 = gesages[x]; 

        keep patid pregdt0379 gesage0379; 
    run;     
%end; 

%mend evw0379; 

%evw0379; 

/*Keep obs where participant is currently pregnant and lmptdt is not missing*/
proc sort data=permhope.evw0380 out=evw0380_sort(keep=patid lmpdt pregnow where=(pregnow=1 and lmpdt ne .));
  by patid; 
run; 

%macro evw0380(); 
/*For currently pregnant patids, keeping pregnancy with lmpdt before entryvisitdt; or before randdt if entryvisit did not happen yet */
data evw0380(drop=entrystratcat2 stratcat2 entryvisdt randdt pregnow); 
    merge evw0380_sort dervhope.master(keep=patid entrystratcat2 stratcat2 entryvisdt randdt);                  
    by patid;                                                                                                                                     
    %if entrystratcat2 ne . %then %do;
      if entrystratcat2=2 and lmpdt <= entryvisdt; 
    %end;
    %else %do;
      if stratcat2=2 and lmpdt <= randdt;
    %end; 
run; 
%mend evw0380; 

%evw0380; 

%macro dates(); 
/*Calculate missing lmpdts and cutoffdts for different stratcats*/
data dates miss_cutoffdt;
    merge dervhope.master(in=inmaster keep=patid entrystratcat2 stratcat2 entryvisdt randdt)
          evw0380(in=in0380)
          evw0379(in=in0379)
          QLW0360(in=in0360)
          anstab(in=inanstab);
    by patid;     
    if inmaster;
    length cutoffdt 5;  

    /*This block of code derives the cutoffdt using the participant information at the time of entryvisit*/
    if entrystratcat2 ne . then do;
        /*If currently pregnant and missing lmpdt and non-missing gesage then do calculation below for lmpdt*/
        if entrystratcat2 in (2) and lmpdt=. and anstab_gesage ne . then lmpdt=entryvisdt-(anstab_gesage*7); 
         /*For the lmpdt calculation, there's a pregdt hierarchy: pregdt0379--> pregdt0360 --> anstab_pregdt*/ 
        else if entrystratcat2 in (3,4) and lmpdt=. then do; 
            if pregdt0379 ne . and gesage0379 ne . then lmpdt=pregdt0379-(gesage0379*7);
                else if pregdt0360 ne . and gesage0360 ne . then lmpdt=pregdt0360-(gesage0360*7);
                else if anstab_pregdt ne . and anstab_gesage ne . then lmpdt=anstab_pregdt-(anstab_gesage*7);
        end; 
        /*If lmpdt missing bc no valid gesage, then assume 40 weeks according to Jessica document and combine with first-nonmissing pregdt*/  
        if lmpdt=. and (pregdt0379 ne . or pregdt0360 ne . or anstab_pregdt ne .) then do; 
            temp_pregdt=coalesce(pregdt0379, pregdt0360, anstab_pregdt); /*Use first non-missing date in this order*/
            if entrystratcat2 in (3,4) then lmpdt=temp_pregdt-(40*7); /*JESSICA document says to use 40 for gesage if no other value is available*/
        end; 
        
        /*Cutoffdt for stratcats 2,3,4 is lmpdt - 6 months*/
        if entrystratcat2 in (2,3,4) and lmpdt ne . then cutoffdt=intnx('month', lmpdt, -6, 'same');     
        /*Cutoffdt for stratcats 1,5 is enrtyvisdt - 1 year*/
        else if entrystratcat2 not in (2,3,4) and entryvisdt ne . then cutoffdt=intnx('year', entryvisdt, -1, 'same');
    end; 

    /*If the cutoffdt is not able to be calculated using information at the time of entryvisit, then the information at the time of enrollment is used*/
    if cutoffdt = . then do;
        /*If currently pregnant and missing lmpdt and non-missing gesage then do calculation below for lmpdt*/
        if stratcat2 in (2) and lmpdt=. and anstab_gesage ne . then lmpdt=randdt-(anstab_gesage*7); 
         /*For the lmpdt calculation, there's a pregdt hierarchy: pregdt0379--> pregdt0360 --> anstab_pregdt*/ 
        else if stratcat2 in (3,4) and lmpdt=. then do; 
            if pregdt0379 ne . and gesage0379 ne . then lmpdt=pregdt0379-(gesage0379*7);
                else if pregdt0360 ne . and gesage0360 ne . then lmpdt=pregdt0360-(gesage0360*7);
                else if anstab_pregdt ne . and anstab_gesage ne . then lmpdt=anstab_pregdt-(anstab_gesage*7);
        end; 
        /*If lmpdt missing bc no valid gesage, then assume 40 weeks according to Jessica document and combine with first-nonmissing pregdt*/  
        if lmpdt=. and (pregdt0379 ne . or pregdt0360 ne . or anstab_pregdt ne .) then do; 
            temp_pregdt=coalesce(pregdt0379, pregdt0360, anstab_pregdt); /*Use first non-missing date in this order*/
            if stratcat2 in (3,4) then lmpdt=temp_pregdt-(40*7); /*JESSICA document says to use 40 for gesage if no other value is available*/
        end; 
        /*Cutoffdt for stratcats 2,3,4 is lmpdt - 6 months*/
        if stratcat2 in (2,3,4) and lmpdt ne . then cutoffdt=intnx('month', lmpdt, -6, 'same');     
        /*Cutoffdt for stratcats 1,5 is enrtyvisdt - 1 year*/
        else if stratcat2 not in (2,3,4) and randdt ne . then cutoffdt=intnx('year', randdt, -1, 'same');
        
        /*This will be temporary, but this flag identifies when a participant had an entryvisit but the info could not be used to calculate cutoffdt so enrollment 
          info was used instead*/
        if entryvisdt ne . and cutoffdt ne . then entry_happen_used_enroll=1;

    end; 

    /*Print patids with no cutoffdts to log and output to separate dataset so they can be investigated further. Print data with cutoffdts
    to dates dataset*/
    if cutoffdt = . then do; 
      output miss_cutoffdt; 
       put 'WARNING: No cutoffdt available for ' 
            patid= entrystratcat2= stratcat2= randdt= entryvisdt= lmpdt= anstab_pregdt= anstab_gesage=
            pregdt0379= pregdt0360= gesage0379= gesage0360=;
    end; 
    else if cutoffdt ne . then output dates;  
    
    label 
      temp_pregdt='First non-missing pregdt (selected in the following order: pregdt0379 pregdt0360 anstab_pregdt) value used as a last resort in lmpdt calc'
      cutoffdt = 'Reflects the date 6 months before lmpdt for patids with stratcat2 = (2,3,4) or reflects the date 1 year before entryvisdt for 
                  patids not in stratcat2 = (2,3,4)';

    format cutoffdt anstab_pregdt pregdt0379 pregdt0360 temp_pregdt date9. ; 
run; 

/*Store number of obs with missing cutoffdts into macro var*/
proc sql noprint; 
  select count(*) into: miss_cutoffdt_obs from miss_cutoffdt; 
quit; 

/*Print obs with missing cutoffdts to CSV for easier readability*/
%if &miss_cutoffdt_obs>0 %then %do; 
  ods csv file="&outrtf./miss_cutoffdt.csv"; 
    proc print data=miss_cutoffdt; 
      var patid entrystratcat2 stratcat2 randdt entryvisdt lmpdt anstab_pregdt anstab_gesage 
        pregdt0379 pregdt0360 gesage0379 gesage0360;
      format anstab_pregdt pregdt0379 pregdt0360 date9. ; 
    run; 
  ods csv close; 
%end; 

%mend dates; 

%dates(); 

/********************************************************************************************************************************************
                                                              DUPLICATES PROGRAMMING                                                  
********************************************************************************************************************************************/

/*Concatenate viral load datasets from the 'AGGREGATING AND STANDARDIZING CD4 DATA FROM ACROSS PROTOCOLS' section to create 
  one comphrensive data view of all viral load information*/
data all_rna; 
    set hope_rna  (in=inhope) smartt_rna (in=insmartt) 
        ampup_rna (in=inamp)  aul_rna    (in=inaul); 
    by patid specdt;
    format quantcd stan_qnt.; 
run;


/*Subset dataset to keep information relevant to HOPE participants and to keep only certain viral load measurements*/
data all_rna_inHOPE(drop=cutoffdt); 
    merge all_rna(in=inall)
          master(in=inmaster)
          dates(in=indates keep=patid cutoffdt)
          copyml_allrslts_flags(in=inflags); 
    by patid; 
    /*Keep patids found in HOPE Master or records where specdts are on or after the cutoffdt*/
    if inmaster or if specdt >= cutoffdt; 
    /*Some patids in HOPE master are not in HOPE viral load CRF so their source value shows up as missing - this corrects that*/
    if source='' then source='HOPE';      
run;


/*Store all duplicates in a dataset that can be looked at to get better understanding of duplicates*/
proc sql;
    create table dups_rna as 
        (select a.*, 
               case 
                   when count(distinct copyml) = 1 then 'All Same'
                   else 'Different'
               end as same_RNA_flag, /*For patid-specdt combos with more than 1 viral load measure, this flag identifies if those measures are all the same or are different*/
               case 
                   when count(*) > 1 then 'Duplicate' 
                   else 'Not Duplicate' 
               end as dup_flag /*This flag identifies if a patid-specdt combo has more than 1 viral load measure*/
        from all_rna_inHOPE as a
        group by patid, specdt
        having count(*) > 1);   
    
    /*Store max number of viral load measures per patid per specdt in a macro variable to be used later*/
    select max(count) into :max_copyml_values
    from (
        select patid, specdt, count(*) as count
        from all_rna_inHOPE
        group by 1,2);
quit;

%macro mkvars(counter);
if counter=&COUNTER then do; 
    copyml_&COUNTER=copyml; quantcd_&COUNTER=quantcd; source_&COUNTER=strip(source); assaysp_&COUNTER=strip(assaysp); assaytyp_&COUNTER=assaytyp; 
end; 
%mend mkvars;

%macro rnaprocess; 

/*Putting all copyml and related data for a given patid specdt combo on 1 line using mkvars macro above*/
  
data rna_onerow; 
    set all_rna_inHOPE(keep=patid specdt copyml quantcd source assaysp assaytyp); 
    by patid specdt;

    /*For %mkvars macro to work, we need to define vars found in macro, retain those vars, and put them into arrays for do-loop processing*/
    length source_1-source_%sysfunc(compress(&max_copyml_values)) $200 assaysp_1-assaysp_%sysfunc(compress(&max_copyml_values)) $70 assaytyp_1-assaytyp_%sysfunc(compress(&max_copyml_values)) 3;  
    retain copyml_1-copyml_%sysfunc(compress(&max_copyml_values)) quantcd_1-quantcd_%sysfunc(compress(&max_copyml_values)) source_1-source_%sysfunc(compress(&max_copyml_values)) counter;
    array retarr_copyml copyml_1-copyml_%sysfunc(compress(&max_copyml_values));
    array retarr_quantcd quantcd_1-quantcd_%sysfunc(compress(&max_copyml_values));
    array retarr_source source_1-source_%sysfunc(compress(&max_copyml_values));
    array retarr_assaysp assaysp_1-assaysp_%sysfunc(compress(&max_copyml_values));
    array retarr_assaytyp assaytyp_1-assaytyp_%sysfunc(compress(&max_copyml_values)); 

    /*Initialize vars*/
    if first.specdt then do;
        counter=1; 
        do i=1 to dim(retarr_copyml);
            retarr_copyml(i)=.;
            retarr_quantcd(i)=.; 
            retarr_source(i)=''; 
            retarr_assaysp(i)='';
            retarr_assaytyp(i)=.; 
        end;
    end;
    else counter+1; 

    /*Call mkvars macro*/
    %do i=1 %to &max_copyml_values; 
      %mkvars(&i);
    %end; 

    if counter>&max_copyml_values then put "ERROR: counter > &max_copyml_values !"  patid= specdt=; 
    
    if last.specdt; /*%mkvars works in such a way that all the related relevant viral load information for a give patid and specdt gets stored on the record where last.specdt is true */

run; 

/*Before viral load measurements can be evaluated, we need to create various flags that will help with the processing*/
data flags_for_process; 

  set rna_onerow; 
  by patid specdt; 
  
  if specdt=. then missdt=1; else missdt=0; /*Flag to identify missing specdts*/
  
  miss_copyml_flag=1; /*Assume a given patid specdt combo has a missing value for all copyml vars*/

  /*Define arrays for do-loop processing*/
  array copymls{*} copyml_1-copyml_%sysfunc(compress(&max_copyml_values)); 
  array quantcds{*} quantcd_1-quantcd_%sysfunc(compress(&max_copyml_values));
  array sources{*} source_1-source_%sysfunc(compress(&max_copyml_values));

  do i=1 to dim(copymls); 
    if copymls[i] ne . then do; 
      miss_copyml_flag=0; /*Set flag to 0 if there is atleast 1 non-missing copyml value for a given patid specdt combo*/
      leave; 
    end; 
  end;

  /* Logic to count number of 'HOPE' copyml values */
  count_hope = 0;
  do i = 1 to dim(sources);
    if sources(i) = 'HOPE' then count_hope + 1;
  end;

  /* Initialize variables and check if all non-missing quantcds are the same */
  all_same_quantcd = 1; /* Assume all are the same until proven otherwise */
  previous_quantcd = .; /* Initialize with a missing value */

  do i = 1 to dim(quantcds);
      if not missing(quantcds[i]) then do;
          if missing(previous_quantcd) then do;
              previous_quantcd = quantcds[i]; /* Set the first non-missing value */
          end;
          else if quantcds[i] ne previous_quantcd then do;
              all_same_quantcd = 0; /* Set flag to 0 if a different value is found */
              leave; /* Exit the loop early */
          end;
      end;
  end;


  /* Initialize flag that identifies if all HOPE quantcds are the same*/
  hope_same_quantcd = 1;
  
  /* Compare HOPE quantcd values */
  do i = 1 to dim(sources);
      if sources[i] = 'HOPE' then do;
          do j = i + 1 to dim(sources);
              if sources[j] = 'HOPE' then do;
                  if quantcds[i] ne quantcds[j] then do;
                      hope_same_quantcd = 0; /* Set flag to 0 if all HOPE quantcds are not the same */
                      leave;
                  end;
              end;
          end;
          if hope_same_quantcd = 0 then leave;
      end;
  end;

/* Initialize flag that identifies if atleast 1 pair of matching quantcds exists. Assume this is not true */
  one_quantcd_match = 0;

  /* Loop through each quantcd and compare with subsequent quantcds */
  do i = 1 to dim(quantcds);
    if not missing(quantcds[i]) then do;
      do j = i + 1 to dim(quantcds);
        if quantcds[i] = quantcds[j] then do;
          one_quantcd_match = i; /*Store the position of the match*/
          leave; /* Exit the inner loop early once a match is found */
        end;
      end;
      if one_quantcd_match then leave; /* Exit the outer loop early if a match is found */
    end;
  end;

  /* Initialize variables */
  previous_copyml = .; /* Initialize with a missing value */
  all_same_copyml = 1; /* Assume all are the same until proven otherwise */

  do i = 1 to dim(copymls);
      if not missing(copymls[i]) then do;
          if missing(previous_copyml) then do;
              previous_copyml = copymls[i]; /* Set the first non-missing value */
          end;
          else if copymls[i] ne previous_copyml then do;
              all_same_copyml = 0; /* Set flag to 0 if different value found */
              leave; /* Exit the loop early if values are not the same */
          end;
      end;
  end;


  /* Initialize flag that identifies if at least 2 HOPE copyml are the same. Assume this is not true */
  hope_same_copyml = 1;
  
  /* Compare HOPE quantcd values */
  do i = 1 to dim(sources);
      if sources[i] = 'HOPE' then do;
          do j = i + 1 to dim(sources);
              if sources[j] = 'HOPE' then do;
                  if copymls[i] ne copymls[j] then do;
                      hope_same_copyml= 0; /* Set flag to 0 if HOPE copymls differ */
                      leave;
                  end;
              end;
          end;
          if hope_same_copyml = 0 then leave;
      end;
  end;    

  /*Initialize flag that identifies if atleast 1 pair of matching copyml exists. Assume this is not true. */
  one_copyml_match = 0;

  /* Loop through each quantcd and compare with subsequent copymls */
  do i = 1 to dim(copymls);
    if not missing(copymls[i]) then do;
      do j = i + 1 to dim(copymls);
        if copymls[i] = copymls[j] then do;
          one_copyml_match = i; /*Store the position of the match*/
          leave; /* Exit the inner loop early once a match is found */
        end;
      end;
      if one_copyml_match then leave; /* Exit the outer loop early if a match is found */
    end;
  end;
  
  /*Flag to identify position of protocols that have matching copymls and quantcds across atleast 2 studies*/
  copyml_quantcd_match = 0; 

  /* Loop through the arrays to check for matching values */
  do i = 1 to dim(copymls);
      do j = i + 1 to dim(copymls);
          if not missing(copymls[i]) and not missing(copymls[j]) then do;
              if copymls[i] = copymls[j] and quantcds[i] = quantcds[j] then do;
                  copyml_quantcd_match = i;
                  leave; /* Exit the inner loop early once a match is found */
              end;
          end;
      end;
  end;  

  /* Declare the hope_pos array to store positions of 'HOPE' */
  array hope_pos[&max_copyml_values] _temporary_;
    
  /* Initialize a counter for the hope_pos array */
  count_hope_pos = 0;

  /* Identify positions of 'HOPE' and count occurrences */
  do i = 1 to dim(sources);
      if sources[i] = 'HOPE' then do;
          count_hope_pos + 1;
          hope_pos[count_hope_pos] = i; /* Store the position of "HOPE" */
      end;
  end;

  /* Store HOPE copyml and quantcd values only*/
  array hope_copymls{&max_copyml_values} _temporary_;
  array hope_quantcds{&max_copyml_values} _temporary_;

  /* Populate arrays with HOPE values */
  do i = 1 to count_hope_pos;
      hope_copymls[i] = copymls[hope_pos[i]];
      hope_quantcds[i] = quantcds[hope_pos[i]];
  end;

  length quantcd_combo $200 symbol $7; /* Initialize quantcd_combo */

  quantcd_combo = ''; /* Initialize quantcd_combo to an empty string */

  do i = 1 to dim(copymls);

    if copymls[i] > 0 then do;
      /* Determine symbol based on quantcd */
      if quantcds[i] = 1 then symbol = '=';
      else if quantcds[i] = 2 then symbol = '>';
      else if quantcds[i] = 3 then symbol = '<';
      else if quantcds[i] = 4 then symbol = 'Unknown';
      else if quantcds[i] = 1.5 then symbol = 'Mean of 2 values with different quantifiers (</=)';
      else if quantcds[i] = . then symbol = ' ';
      else symbol = '/';

      /* Add symbol to quantcd_combo with '/' separator */
      if quantcd_combo ne '' then quantcd_combo = catx('/', quantcd_combo, strip(symbol));
      else quantcd_combo = strip(symbol);

    end;
  end;

  /* Remove leading and trailing spaces */
  quantcd_combo = strip(quantcd_combo);

  label missdt = 'Flag for Missing specdt'
        miss_copyml_flag = 'Flag for obs with all Missing copyml'
        count_hope = 'Count of HOPE values in an observation'
        all_same_quantcd = 'Flag for if non-missing quantcds for a given obs are all the same'
        previous_quantcd = 'Previous quantcd value for comparison'
        hope_same_quantcd = 'Flag to identify if atleast 2 HOPE quantcds are the same'
        one_quantcd_match = 'Flag for at least one pair of matching quantcds'
        previous_copyml = 'Previous copyml value for comparison'
        all_same_copyml = 'Flag for if all copyml for an observation are the same'
        hope_same_copyml = 'Flag to identify if atleast 2 HOPE copymls are the same'
        one_copyml_match = 'Flag for at least one pair of matching copymls'
        copyml_quantcd_match = 'Flag for matching copyml and quantcd values'
        count_hope_pos = 'Count of positions of HOPE in sources'
        quantcd_combo = 'Combined quantcd symbols'
        symbol = 'Symbol for quantcd';
run; 

data keepone; 
  
  set flags_for_process; 
  by patid specdt; 
  length copyml 5 quantcd 3 source $200 dataqcat 8 dataqsp $200.;

  /*Define arrays for do-loop processing*/
  array copymls{*} copyml_1-copyml_%sysfunc(compress(&max_copyml_values)); 
  array quantcds{*} quantcd_1-quantcd_%sysfunc(compress(&max_copyml_values));
  array sources{*} source_1-source_%sysfunc(compress(&max_copyml_values));

  /*If counter=1 then these arent duplicates and minimal processing is required*/
  if counter=1 then do; 
    copyml=coalesce(of copymls[*]); 
    quantcd=quantcds[whichn(copyml, of copymls[*])];
    source=sources[whichn(copyml, of copymls[*])]; 
    dataqcat=1; dataqsp='No duplicate: use all values as are'; 
  end; 
  /*Begin processing duplicates*/
  else if counter>1 and count_hope = 1 then do; 
    *********************************************************************************;    
    **If HOPE and non-HOPE defer to HOPE values**
    *********************************************************************************;
      copyml=copymls[whichc('HOPE', of sources[*])];
      quantcd=quantcds[whichc('HOPE', of sources[*])];
      dataqcat=1; dataqsp='HOPE and Non-HOPE values: Use HOPE values'; source = 'HOPE'; 
      *********************************************************************************;    
      ** if more than 1 viral load measure and one of the measures is from HOPE and the 
         discrepancy is > 10, then still use HOPE values but dataqcat is 'Fixed' **
      *********************************************************************************;
      if range(of copyml_1-copyml_%sysfunc(compress(&max_copyml_values))) > 10 then dataqcat=2;

  end; 
  else if count_hope > 1 then do; 

    source = 'HOPE'; 
    dataqcat = 2; 

    /* Declare the hope_pos array to store positions of 'HOPE' */
    array hope_pos[&max_copyml_values] _temporary_;
    
    /* Initialize a counter for the hope_pos array */
    count_hope_pos = 0;

    /* Identify positions of 'HOPE' and count occurrences */
    do i = 1 to dim(sources);
        if sources[i] = 'HOPE' then do;
            count_hope_pos + 1;
            hope_pos[count_hope_pos] = i; /* Store the position of 'HOPE' */
        end;
    end;

    match_found = 0; /* Flag to identify if one of the HOPE values matches a non-HOPE value */
    do i = 1 to count_hope_pos;
        hope_copyml = copymls[hope_pos[i]];  /* Use hope_pos array to ensure the copyml value associated with HOPE from the copymls array is assigned to hope_copyml */   
        do j = 1 to dim(sources);
            if sources[j] ne 'HOPE' and copymls[j] = hope_copyml then do;
                match_found = hope_pos[i]; /* Store the position of the matching HOPE value */
                leave; /* Exit the loop once a match is found */
            end;
        end;
        if match_found then leave;
    end;

    /* Store HOPE copyml values only in order to calculate min, max, and geometric mean for HOPE values only */
    array hope_copymls{&max_copyml_values} _temporary_;
    array hope_quantcds{&max_copyml_values} _temporary_;

    /* Populate arrays with HOPE values */
    do i = 1 to count_hope_pos;
        hope_copymls[i] = copymls[hope_pos[i]];
        hope_quantcds[i] = quantcds[hope_pos[i]];
    end; 

    if match_found=0 then do; 
    /*If match found between HOPE copyml and non-HOPE copyml*/
      if hope_same_quantcd = 1 and hope_same_copyml = 0 then do;
      /*If all the HOPE quantcds are the same and the HOPE copyml differ */  
        quantcd = coalesce(of hope_quantcds[*]);

        /* Determine copyml based on quantcd */
        if quantcd = 3 then do;
          copyml = min(of hope_copymls[*]); dataqsp = 'RNA differs, quantcd < - use lower'; 
        end; 
        else if quantcd = 2 then do ;
          copyml = max(of hope_copymls[*]); dataqsp = 'RNA differs, quantcd > - use higher'; 
        end; 
        else if quantcd = 1 then do; 
        /*log_sum,non_zero_count,on_zero_value are all vars used for geomean calculation*/
          log_sum = 0;                      
          non_zero_count = 0; 
          non_zero_value = .;
          do i = 1 to dim(hope_copymls);
              if hope_copymls[i] > 0 then do;
                log_sum + log(hope_copymls[i]);
                non_zero_count + 1;
                non_zero_value = hope_copymls[i];
              end;
          end;
          if non_zero_count >= 2 then do;
            copyml = exp(log_sum / non_zero_count);
            dataqsp = cats('RNA differ, quantcd =, - use geomean');
          end;
          else if non_zero_count = 1 then do;
            copyml = non_zero_value;
            dataqsp = 'Single non-zero RNA value from HOPE - use value';
          end;
        end; 
      end; 
      else if hope_same_quantcd = 0 and hope_same_copyml = 0 then do;
      /*If all the HOPE quantcds are different and the HOPE copyml differ */  
        log_sum = 0;
        non_zero_count = 0;
        non_zero_value = .;
        do i = 1 to dim(hope_copymls);
          if hope_copymls[i] > 0 then do;
              log_sum + log(hope_copymls[i]);
              non_zero_count + 1;
              non_zero_value = hope_copymls[i];
          end;
        end;

        if non_zero_count >= 2 then do;
          copyml = exp(log_sum / non_zero_count);
          if countw(cats(of hope_quantcds[*]), '1') > 0 and countw(cats(of hope_quantcds[*]), '3') > 0 then do;  
          *****************************************************************************************;
          ** If </= or =/< found then use </= **;
          *****************************************************************************************;
              quantcd = 1.5; dataqsp = cats('RNA differ, quantcds (', quantcd_combo, ') - take geomean and use </='); 
          end; 
          else if countw(cats(of hope_quantcds[*]), '1') > 0 then do; 
          *****************************************************************************************;
          ** If = is found with > or ' ' then use = **;
          *****************************************************************************************;
              quantcd = 1; dataqsp = cats('RNA differ, quantcds (', quantcd_combo, ') - take geomean and use ='); 
          end;
          else if countw(cats(of hope_quantcds[*]), '3') > 0 then do; 
          *****************************************************************************************;
          ** If < is found with > or ' ' then use < **;
          *****************************************************************************************;
              quantcd = 3; dataqsp = cats('RNA differ, quantcds (', quantcd_combo, ') - take geomean and use <'); 
          end;
          else if countw(cats(of hope_quantcds[*]), '2') > 0 then do; 
          *****************************************************************************************;
          ** If  /> or >/  found then use > **;
          *****************************************************************************************;
            quantcd = 2; dataqsp = cats('RNA differ, quantcds (', quantcd_combo, ') - take geomean and use >'); 
          end;
        end;
        else if non_zero_count = 1 then do;
          copyml = non_zero_value;
          quantcd = hope_quantcds[whichn(non_zero_value, of hope_copymls[*])];
          dataqsp = 'Single non-zero RNA value from HOPE - use value';
        end;
      end;  
      else if hope_same_quantcd = 0 and hope_same_copyml = 1 then do;
      *****************************************************************************************;
      ** If viral loads are the same but quantcds differ **;
      *****************************************************************************************;  
        copyml = coalesce(of hope_copymls[*]); 
        if countw(cats(of hope_quantcds[*]), '1') > 0 and countw(cats(of hope_quantcds[*]), '3') > 0 and copyml not in (20, 40) then do; 
        *************************************************************************************************;
        ** If </= or =/< found and copyml not = to common assay cutoff values of 20 or 40, then use = **;
        *************************************************************************************************;
            quantcd = 1; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use ='); 
        end; 
        else if countw(cats(of hope_quantcds[*]), '1') > 0 and countw(cats(of hope_quantcds[*]), '3') > 0 and copyml in (20, 40) then do; 
        *****************************************************************************************;
        ** If </= or =/< found and copyml = a common assay test cutoff then use <**;
        *****************************************************************************************;
            quantcd = 3; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use <'); 
        end; 
        else if countw(cats(of hope_quantcds[*]), '1') > 0 then do; 
        *****************************************************************************************;
        ** If =/> or =/' ' found then use = **;
        *****************************************************************************************;
            quantcd = 1; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use ='); 
        end;
        else if countw(cats(of hope_quantcds[*]), '3') > 0 then do; 
        *****************************************************************************************;
        ** If  >/< or </' '  found then use < **;
        *****************************************************************************************;
            quantcd = 3; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use <'); 
        end;
        else if countw(cats(of hope_quantcds[*]), '2') > 0 then do; 
        *****************************************************************************************;
        ** If  > found with missing or unknown then use > **;
        *****************************************************************************************;
            quantcd = 2; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use >'); 
        end;
      end; 
    end;
    else if match_found>0 then do; 
      copyml=copymls[hope_pos[match_found]];  quantcd=quantcds[hope_pos[match_found]];
      dataqsp = 'Multiple HOPE values; 1 matches a non-HOPE value - Defer to matching HOPE value'; 
    end; 
  end; 
  else if counter > 1 and count_hope = 0 then do; 
  ***********************************************;
  ** If non HOPE duplicates **;
  ***********************************************;

    dataqcat = 2; /* Set data quality to fixed */

    if copyml_quantcd_match > 0 then do; 
    **********************************************************;
    ** If copyml and quantcd from one study matches copyml and 
       quantcd from another study - use values from either study**;
    **********************************************************;
      copyml = copymls[copyml_quantcd_match]; quantcd = quantcds[copyml_quantcd_match]; 
      source = sources[copyml_quantcd_match]; dataqsp='RNA and quantcd from one study match the set from another study - defer to either set';
    end;
    else if all_same_quantcd=1 and all_same_copyml = 0 then do; 
    **********************************************************;
    ** If all the quantcds are the same and viral load differ **;
    **********************************************************;
        
      quantcd = coalesce(of quantcds[*]);
      
      *****************************************************************************************;
      ** Take lower when quantcd is <, greater when it is > and geomean when it is = and two or 
         more non-zero values**;
      *****************************************************************************************;
        if quantcd = 3 then do;
            copyml = min(of copymls[*]);
            dataqsp = 'RNA differs, qnt < - use lower'; 
            source = sources[whichn(copyml, of copymls[*])];
        end; 
        else if quantcd = 2 then do ;
            copyml = max(of copymls[*]);
            dataqsp = 'RNA differs, qnt > - use higher'; 
            source = sources[whichn(copyml, of copymls[*])];
        end; 
        else if quantcd = 1 then do; 
            log_sum = 0;
            non_zero_count = 0;
            non_zero_value = .;
            do i = 1 to dim(copymls);
                if copymls[i] > 0 then do;
                    log_sum + log(copymls[i]);
                    non_zero_count + 1;
                    non_zero_value = copymls[i];
                end;
            end;

            if non_zero_count >= 2 then do;
                copyml = exp(log_sum / non_zero_count);
                dataqsp = cats('RNA differ, quantcd = - use geomean');
                if source = '' then source =strip(sources[i]);
                else source = catx('/ ', strip(source), sources[i]);
            end;
        end; 
        else if non_zero_count = 1 then do; 
            copyml = non_zero_value;
            dataqsp = 'Single non-zero RNA value - use value';
            source = sources[whichn(copyml, of copymls[*])];
        end; 
      end;
    end;  
    else if all_same_quantcd=0 and all_same_copyml = 1 then do; 
    *****************************************************************************************;
    ** If all the quantcds are not the same and all viral load measures are the same **;
    *****************************************************************************************;
      copyml=coalesce(of copymls[*]); 

      if countw(cats(of quantcds[*]), '1') > 0 and countw(cats(of quantcds[*]), '3') > 0 and copyml not in (20,40) then do; 
      *****************************************************************************************;
      ** If < and = found and copyml ne to a common assay cutoff then use = **;
      *****************************************************************************************;
          quantcd = 1; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use =');
      end; 
      else if countw(cats(of quantcds[*]), '1') > 0 and countw(cats(of quantcds[*]), '3') > 0 and copyml in (20,40) then do; 
      *****************************************************************************************;
      ** If < and = found and copyml eq to a common assay cutoff then use < **;
      *****************************************************************************************;
          quantcd = 3; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use <');
      end; 
      else if countw(cats(of quantcds[*]), '1') > 0 then do;
      *****************************************************************************************;
      ** If = and > found or = and ' ' found then use = **;
      *****************************************************************************************;    
          quantcd = 1; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use =');
      end; 
      else if countw(cats(of quantcds[*]), '3') > 0 then do;
      *****************************************************************************************;
      ** If </> or ' ' /< found then use < **;
      *****************************************************************************************; 
          quantcd = 3; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use <');
      end; 
      else if countw(cats(of quantcds[*]), '2') > 0 then do;
      *****************************************************************************************;
      ** If > /' ' or ' ' /> found then use > **;
      *****************************************************************************************;   
          quantcd = 2; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use >');
      end; 

      source = sources[whichn(copyml, of copymls[*])];

    end;    
    else if all_same_quantcd=0 and all_same_copyml = 0 then do; 
    *****************************************************************************************;
    ** If all the quantcds are not the same and all viral load measures are not the same **;
    *****************************************************************************************;
        dataqcat = 2;   
        log_sum = 0;
        non_zero_count = 0;
        non_zero_value = .;
        do i = 1 to dim(copymls);
            if copymls[i] > 0 then do;
                log_sum + log(copymls[i]);
                non_zero_count + 1;
                non_zero_value = copymls[i];
                if source ne '' then source = catx('/', strip(source), sources[i]);
                else source = strip(sources[i]);
            end;
        end;

        if non_zero_count >= 2 then do;
            copyml = exp(log_sum / non_zero_count);
            if countw(cats(of hope_quantcds[*]), '1') > 0 and countw(cats(of hope_quantcds[*]), '3') > 0 then do;  
            *****************************************************************************************;
            ** If </= or =/< found then use </= **;
            *****************************************************************************************;
                quantcd = 1.5; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use geomean'); 
            end; 
            else if countw(cats(of hope_quantcds[*]), '1') > 0 then do; 
            *****************************************************************************************;
            ** If = is found with > or '' then use = **;
            *****************************************************************************************;
                quantcd = 1; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - take geomean and use ='); 
            end;
            else if countw(cats(of hope_quantcds[*]), '3') > 0 then do; 
            *****************************************************************************************;
            ** If  >/< or </'' found then use < **;
            *****************************************************************************************;
                quantcd = 3; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use <'); 
            end;
            else if countw(cats(of hope_quantcds[*]), '2') > 0 then do; 
            *****************************************************************************************;
            ** If ''/> or >/''  found then use > **;
            *****************************************************************************************;
                quantcd = 2; dataqsp = cats('RNA same, quantcds (', quantcd_combo, ') - use >'); 
            end;
        end;
        else if non_zero_count = 1 then do;
            copyml = non_zero_value;
            quantcd = quantcds[whichn(non_zero_value, of copymls[*])];
            dataqsp = 'Single non-zero RNA value - use that value';
            source = sources[whichn(non_zero_value, of copymls[*])];
        end; 
    end; 
  else if missdt=0 and miss_copyml_flag=0 and counter>1 then do; 
    put "WARNING: unexpected scenario!" patid= specdt= counter= source_1= source_2=
                                        quantcd_1=stan_qnt. copyml_1= quantcd_2=stan_qnt. copyml_2= 
                                        quantcd_3=stan_qnt. copyml_3= ;
    
    warn=1; /*We only want a warning when if-thens cannot evaluate a scemario that has all the necessary components: specdt and atleast 1 copyml value*/
  end;

  format dataqcat DATAQUALF. quantcd stan_qnt.;
  label
    dataqsp='Data quality description'
    dataqcat='Data quality indicator';      

run; 

/********************************************************************************************************************************************
                                                   VIRAL LOAD RELATED VARIABLES PROGRAMMING                                                  
********************************************************************************************************************************************/

/*Finalizing assaysp, assaytyp, and logrna variables*/
data keepone; 
    set keepone; 

    length assaysp $70 assaytyp 3;
    /*Define arrays for do-loop processing*/
    array sources {*} $ source_1-source_%sysfunc(compress(&max_copyml_values)); 
    array assaysps {*} $ assaysp_1-assaysp_%sysfunc(compress(&max_copyml_values));
    array assaytyps {*} assaytyp_1-assaytyp_%sysfunc(compress(&max_copyml_values));


    if counter=1 then do;  
        assaysp=coalescec(of assaysps[*]); 
        assaytyp=coalesce(of assaytyps[*]);
    end; 
    else if count_hope > 0 then do; 
    *********************************************************************************;
    **Any record where dataqsp has HOPE in it, it means the HOPE value was used so 
      source, assaysp, and assaytyp need to be assigned the values that were associated
      with HOPE**
    *********************************************************************************;
        assaysp  = assaysps[whichc('HOPE', of sources[*])];
        assaytyp = assaytyps[whichc('HOPE', of sources[*])]; 
    end; 
    else if (find(dataqsp, 'geomean') > 0 and whichc('PH400', of sources[*]) > 0) or source = 'PH400' then do; 
    *********************************************************************************;
    **Aside from HOPE, PH400 is the only study that reports assaysp and assaytyp. If 
      an obs has a viral load value that is the mean of viral load measures from 
      across studies(excluding HOPE), and one of the studies contributing a copyml 
      value was PH400, then we can just use the assaysp and assaytyp from PH400**
    *********************************************************************************;    
        assaysp  = assaysps[whichc('PH400', of sources[*])];
        assaytyp = assaytyps[whichc('PH400', of sources[*])]; 
    end; 
    else do; 
    *********************************************************************************; 
    *If copyml is geomean of copyml from SMARTT and AMP Up, then both studies do not 
        have assaysp or assaytyp values so we have them as blanks*
    *********************************************************************************; 
        assaysp=''; 
        assaytyp=.; 
    end; 
    
  
    /*************************************************************************************
    AMONG THE 3 LOGRNA DERIVATION OPTIONS(AMP UP, HOPE, SMARTT), JESSICA DEB AND KATHY DECIDED 
    SMARTT DERIVATION WAS TO BE USED. AMP UP AND HOPE DERIVATIONS CAN BE FOUND BELOW IN 
    CASE INTERESTED
    
    AMP UP(comb_rnamrg) logrna derivation: if copyml ne 0 then logrna=log10(copyml);
    
    HOPE(master) logrna derivation: 
        if rnadt ne . then do;
            *LLD correction and transform to log rna;
            if quantcd=1 and copyml>0 then copyml_lld=copyml-1;
                else copyml_lld=copyml;length source $200; 
    source=strip('PH400');

        if copyml_lld>0 then logrna=log10(copyml_lld);
            else if copyml_lld=0 then logrna=0;
    end;
    
  *************************************************************************************/
    length logrna 8; 
    if copyml=0 then logrna=log10(1);
    else logrna=log10(copyml);

run; 

%mend rnaprocess;

%rnaprocess; 

/********************************************************************************************************************************************
                                                  FINAL PROCESSING AND CHECKS                                                 
********************************************************************************************************************************************/

/*Create final dataset with ammended labels*/
data rna;
    set keepone; 
    label 
       assaysp  = 'Specify other assay code'
       assaytyp = 'Assay type code'
       copyml   = 'HIV RNA (copies/mL) - (geo mean if > 1 value per specdt for non-HOPE measures)'
       logrna   = 'log10 HIV RNA (copies/mL) - (0 if copyml=0)'
       quantcd  = 'Quantifier'
       source   = 'Protocols copyml value derived from'
       specdt   = 'Specimen date' ;

    keep patid specdt source logrna copyml quantcd dataqsp dataqcat assaysp assaytyp; 
run; 


/********************************************************************************************************************************
                                    CREATING INDICATOR VARS PERTAINING TO RNA DATA COLLECTION 
********************************************************************************************************************************/

/* Sort the dataset hxw0233 by patid and visit date, keeping only necessary variables: patid, visitdt, and norslt indicators */
proc sort data=permhope.hxw0233 out=hxw0233_sort(keep=patid visitdt allrslts norslt2 norslt3 norslt4); 
    by patid visitdt; 
run; 

/* Create a new dataset "indicators" with renamed indicator variables 
   for clearer interpretation. For each patid, retain only the latest
   visit date (last.patid)*/
data indicators(rename=(norslt2=Pending_SomeRecords 
                        norslt3=Pending_AllRecords 
                        norslt4=Recs_NoLongerExpected)); 
    set hxw0233_sort; 
    by patid visitdt; 
    /*If PID answered Yes to allrslts (which means they skipped the questions that correspond to norslt2-norslt4), 
      then these indicator variables should be set to No instead of left Unknown*/
    if allrslts = 1 then do;
        norslt2 = 2; 
        norslt3 = 2; 
        norslt4 = 2; 
    end; 
    if last.patid;  
    drop allrslts; 
run; 

/* Merge the indicators dataset with the RNA dataset to add the renamed 
   indicators (Pending_SomeRecords, Pending_AllRecords, Recs_Unreceived_NoLongerExpected) 
   based on the patid. Only keep records present in the RNA dataset.*/
data dervhope.rna; 
    merge rna(in=ina) 
          indicators(keep=patid Pending_SomeRecords 
                                Pending_AllRecords 
                                Recs_NoLongerExpected); 
    by patid; 
    if ina; /* Keeps records from the RNA dataset only */
    label Pending_SomeRecords = "Pending some medical records at most recent form entry - Yes/No"
          Pending_AllRecords = "Pending all medical records at most recent form entry - Yes/No"
          Recs_NoLongerExpected = "Outstanding requested medical records not received and no longer expected at most recent form entry - Yes/No";
run; 

%pgmlabel(dervhope.rna, LABEL=%str(Longitudinal Viral Load History for HOPE));

 

  