#!/bin/bash
# Author: Etienne CAMENEN
# Date: 2018
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
# Short description: performs multi-variate analysis (PCA, CCA, PLS, RGCCA) and projects the variables and samples into a bi-dimensional space.

#Settings files
FILE_ERR="false"
OUTFILES=( 'individuals.pdf' 'corcircle.pdf' 'top_variables.pdf' 'ave.pdf' 'design.pdf' 'individuals.tsv' 'variables.tsv' 'rgcca.result.RData')

#Initialization
declare -x INFILE FUNC OPAR WARN
declare -i PARAMETER NBFAIL=0 NBTEST=0 EXIT
declare -a TESTS WARNS EXITS
echo '' > resultRuns.log
echo '' > warnings.log
export LANG=en_US.utf8

setUp(){
    INFILE="data/agriculture.tsv,data/industry.tsv,data/politic.tsv"
    EXIT=0
    PARAMETER=0
    WARN="fgHj4yh"
    WARNS=()
    FUNC=${FUNCNAME[1]}
    TESTS=()
    EXITS=()
    printf "\n- ${FUNC}: "
}

tearDown(){
    rm -r temp/
    mkdir temp/
}

########### ERRORS CATCH ###########

testError(){
    local BOOLEAN_ERR="false"
    local MSG=""
    local ACTUAL_OUTPUT=$(cat temp/log | tr '\n' ' ' )

    [ $1 -ne ${EXIT} ] && {
        MSG=${MSG}"Program exited with bad error code: $1.\n"
        BOOLEAN_ERR="true"
    }

     if [ $1 -eq 0 ]
     then

        for i in ${OUTFILES[@]}; do
	    	testFileExist ${i}
	    done

	    [ ${FILE_ERR} == "true" ] && {
           MSG=${MSG}"not created.\n"
           BOOLEAN_ERR="true"
        }

     else

        if [[ ${ACTUAL_OUTPUT} != *"$WARN"* ]]; then
            MSG=${MSG}"Expected warnings not found.\n"
            echo "Expected: $WARN\n" >> warnings.log
            echo "Actual: $ACTUAL_OUTPUT" >> warnings.log
            BOOLEAN_ERR="true"
        fi
    fi

    [ ${BOOLEAN_ERR} == "true" ] && {
	    ERRORS=${ERRORS}"\n##Test \"${TESTS[$2]}\" in $FUNC: \n$MSG"
	    return 1
    }
    return 0
}

testFileExist(){

   [ ! -f "temp/"$1 ]  && {
        MSG=${MSG}"$1 "
        FILE_ERR="true"
    }
}

printError(){
    testError $@
    if [ $? -ne 0 ]; then
        echo -n "E"
        let NBFAIL+=1
    else echo -n "."
    fi
}

########### RUN PROCESS ###########

run(){
    let NBTEST+=1
	printf "\n\n$NBTEST. ${TESTS[$PARAMETER]}\n" >> resultRuns.log 2>&1
    let PARAMETER+=1
    Rscript inst/launcher.R -d ${INFILE} ${O_PAR} $@ > temp/log 2>&1
}

getElapsedTime(){
    local END_TIME=$(date -u -d $(date +"%H:%M:%S") +"%s")
    local ELAPSED_TIME=$(date -u -d "0 $END_TIME sec - $1 sec" +"%H:%M:%S")
    echo "Time to run the process ${ELAPSED_TIME:0:2}h ${ELAPSED_TIME:3:2}min ${ELAPSED_TIME:6:2}s"
}

setOutputPar(){
    OPAR=""
    for i in `seq 0 $((${#OUTFILES[@]} -1))`; do
        let j=${i}+1
        O_PAR=${O_PAR}"--o"${j}" temp/"${OUTFILES[i]}" "
    done
}

########### TESTS ###########

test(){
    for i in `seq 0 $((${#TESTS[@]} -1))`; do
        run "-d $INFILE ${TESTS[i]}" > temp/log 2>&1
        local ACTUAL_EXIT=$?
        cat temp/log >> resultRuns.log

        if [ ! -z $1 ]; then
            WARN="${WARNS[i]}"
            if [ $1 == "2" ]; then
                EXIT="${EXITS[i]}"
            fi
        fi

        printError ${ACTUAL_EXIT} ${i}

        tearDown
    done
}

testsDefault(){
    setUp
    TESTS=( '' )
    test
}

testsBlocksNrow(){
    cat data/agriculture.tsv | head -n -1 > temp/agriculture2.tsv
    setUp
    INFILE="temp/agriculture2.tsv,data/industry.tsv,data/politic.tsv"
    TESTS=( '' )
    test
}

testsSep(){
    setUp
    TESTS=( '--separator 1')
    test
}

testsSepBad(){
    setUp
    EXIT=1
    WARN="separator should be comprise between 1 and 3 (1: Tabulation, 2: Semicolon, 3: Comma) [by default: 2], not be equal to"
    WARNS=( "${WARN} 0" "${WARN} 4" "agriculture.tsv has an only-column. Check the separator." )
    TESTS=( '--separator 0' '--separator 4' '--separator 2' )
    test 1
}

testsScheme(){
    setUp
    for i in `seq 0 3`; do
        let j=${i}+1
        TESTS[i]='--scheme '${j}
    done
    test
}

testsSchemeBad(){
    setUp
    EXIT=1
    WARN="scheme should be comprise between 1 and 4 [by default: 2], not be equal to"
    WARNS=( "$WARN 0." "$WARN 5.")
    TESTS=( '--scheme 0' '--scheme 5' )
    test 1
}

testsResponse(){
    setUp
    TESTS=( '--group data/political_system.tsv' '--group data/response2.tsv' '--group data/response3.tsv' )
    test
}

testsResponseBad(){
    cat data/political_system.tsv | head -n -1 > temp2/response.tsv
    setUp
    EXITS=(0 1 0)
    WARNS=( "" "test.tsv does not exist" "There is multiple columns in the response file. By default, only the first one is taken in account.")
    TESTS=( '--group temp2/response.tsv' '--group test.tsv' "--group data/agriculture.tsv")
    test 2
}

testsConnection(){
    setUp
    TESTS=( '--superblock -c data/connection.tsv')
    test
}

testsConnectionBad(){
    setUp
    EXIT=1
    WARNS=( "The connection file should contain only 0 or 1." "The connection file should be a symmetric matrix." "The connection file should not contain only 0." "connection matrix should have the same number of columns (actually 4) than the number of blocks (3)." )
    cat data/connection2.tsv | tr '[1]' '[2]' > temp2/connection.tsv
    cat data/connection2.tsv | head -n -1 > temp2/connection2.tsv
    cat data/connection2.tsv | head -n -1 | cut -f -3  > temp2/connection3.tsv
    TESTS=( '--superblock -c temp2/connection.tsv' '--superblock -c temp2/connection2.tsv' '--superblock -c temp2/connection3.tsv' "--superblock -c data/connection2.tsv")
    test 1
}

testHeader(){
    setUp
    WARN="agriculture file contains qualitative data. Please, transform them in a disjunctive table. Possible mistake: header parameter is disabled, check if the file doesn't have one."
    TESTS=( '-H --group data/response3.tsv' )
    test
}

testExcel(){
    setUp
    INFILE="data/blocks.xlsx"
    TESTS=( '' )
    test
}

testFileCharacter(){
    paste data/agriculture.tsv data/political_system.tsv > temp/dummyFile.tsv
    setUp
    EXIT=1
    WARN="dummyFile contains qualitative data. Please, transform them in a disjunctive table."
    INFILE="temp/dummyFile.tsv"
    TESTS=( '' )
    test
}

testTauBad(){
    setUp
    EXIT=1
    WARN="tau should be comprise between 0 and 1 or should correspond to the character 'optimal' for automatic setting (currently equals to"
    WARNS=("$WARN 2).")
    TESTS=('--tau 2')
    test 1
}

testTau(){
    setUp
    j=0
    for i in `seq 0 0.2 1 `; do
        TESTS[${j}]='--tau '${i}
        let j=${j}+1
    done
    TESTS[6]="--tau optimal"
    TESTS[7]="--tau 0,1,0,0.75"
    TESTS[8]="--tau 0.1"
    test
}

testOtherRGCCAparam(){
    setUp
    TESTS=( '--scale' '--superblock' )
    test
}

testNmark(){
    setUp
    TESTS=('--nmark 2' '--nmark 100' '--nmark 1000000')
    test
}

testNmarkBad(){
    setUp
    EXIT=1
    WARN="nmark should be higher than or equal to 2."
    WARNS=( "$WARN" "$WARN" )
    TESTS=('--nmark 1' '--nmark 0')
    test 1
}

testsInit(){
    setUp
    for i in `seq 0 1`; do
        let j=${i}+1
        TESTS[i]='--init '${j}
    done
    test
}

testsInitBad(){
    setUp
    EXIT=1
    WARN='init should be comprise between 1 and 2 (1: Singular Value Decompostion , 2: random) [by default: SVD].'
    TESTS=( '--init 0' '--init 3' '--init 0.3')
    test
}

testNcomp(){
    setUp
    TESTS=( '--ncomp 2' '--ncomp 2,2,2' '--compx 1' '--compy 2')
    test
}

testNcompBad(){
    setUp
    EXIT=1
    WARN='ncomp should be higher than or equal to 2.'
    WARNS=("$WARN" "$WARN" "ncomp should be comprise between 2 and 2, the number of variables of the block (currently equals to 4)" "$WARN" "ncomp should be numeric.")
    TESTS=( '--ncomp 0' '--ncomp 1' '--ncomp 2,4,2,2' '--ncomp 2,0.2,2,2' '--ncomp kjlj')
    test 1
}

testNcompXY(){
    setUp
    TESTS=( '--compx 1' '--compy 2' '--compx 2.1')
    test
}

testNcompXYBad(){
    setUp
    EXIT=1
    WARN="should be higher than or equal to 1."
    WARN0='is currently equals to'
    WARN1='and should be comprise between 1 and'
    WARN2='(the number of component for the selected block).'
    WARNS=("compx $WARN" "compy $WARN" "compx $WARN0 3 $WARN1 2 $WARN2" "compy $WARN0 3 $WARN1 2 $WARN2" 'integer expected, got "mlkmk"')
    TESTS=( '--compx 0' '--compy 0' '--compx 3' '--compy 3' '--compy mlkmk' )
    test 1
}

testBlock(){
    setUp
    for i in `seq 0 3`; do
        TESTS[i]='--block '${i}
    done
    test
}

testBlockBad(){
    setUp
    EXIT=1
    WARN='block should be lower than 4 (the maximum number of blocks), not be equal to 100.'
    TESTS=( '--block 100' )
    test
}

########### MAIN ###########

START_TIME=$(date -u -d $(date +"%H:%M:%S") +"%s")
mkdir temp/ temp2/
setOutputPar

testsDefault
testsSep
testsSepBad
testsScheme
testsSchemeBad
testsResponse
testsResponseBad
testsConnection
testsConnectionBad
testHeader
#testExcel
testFileCharacter
testsBlocksNrow
testTau
testTauBad
testOtherRGCCAparam
testNmark
testNmarkBad
#testsInit
#testsInitBad
testNcomp
testNcompBad
testNcompXY
testNcompXYBad
testBlock
testBlockBad

rm -r temp/ temp2/
printf "\n$NBTEST tests, $NBFAIL failed.$ERRORS\n"
getElapsedTime ${START_TIME}
[[ -z $(cat warnings.log) ]] && rm warnings.log
[[ -z ${ERRORS} ]] || exit 1
exit 0