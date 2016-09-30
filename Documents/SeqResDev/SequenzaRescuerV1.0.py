# This script is used to rescue sequenca output
# This is a WIP

# Filter name: "/home/rdmorin/agilent_sureselect_all_exons_v5_and_utr.sort.merge.bed"
import sys
import subprocess
import numpy
import os

# Obtains arguments from the command line, stores them in an array, and returns the array.
# If there are no arguments, prints out a list of all availible options
def GetStartingVariables():
	"Obtains parameters from the command line"


	print("")
	# If there are no arguments, print out the help message
	if len(sys.argv) == 1 :
	
		print("This script is designed to rescue variants that would normally be missed in the sequqnza cleanup")
		print("The following commands are availible\n")
		print("     -f <filterSeq.fa>    Reference file used in mpileup")
		print("     -d <int>             Minimum filter depth [Default:8]")
		print("     -D <int>             Maximum filter depth PERCENTILE [Default:95]")
		print("     -l <filterFile.bed>  Only look for variants sites")
		print("     -o <amendment>       Ammends file name with this string [Default:_filt]")
		print("     -p <path>            Specify an output path for the files [Default is inputpath]")
		print("     -debug               Prints out all the variables used for mpileup,call, and filter, as well as input and outut file names")
		print("     -q <int>             Minimum quality read percentage [Default:5]")
		print("")
		exit()
	else:
		cmdArgs=[]

		for i in range(1, len(sys.argv)):

			cmdArgs.append(sys.argv[i])

	return cmdArgs

# Prints the output file names and parameters listed in the dictionary to the console
def printRawParameters(parameters):

	print("\nOutputFileNames")
	for names in parameters["output"]:
		print(names)


	print("\nDictionary parameters")
	for key, value in parameters.items():
		print(key + ": " + str(value))

def setOutputFileNames(inputSequences, ammendedPath, ammendedName="_filt"):
	"This sets an output file name for an input file using defaults or user variables (if supplied)"

	# If no input files were provided, quit the program
	outputNames=[]
	if len(inputSequences) == 0:
		print("No input file was specified")
		exit()

	# If there are input files,
	else:

		fileNumber=1
		# For each file provided
		for inputName in inputSequences:
			
			# Obtains the start of the file name in the string (Removed directory bit)
			inputFileStart=len(inputName)
			for char in reversed(inputName):

				if char == "/":

					break
				else:

					inputFileStart-=1

			# Obtains the end of the file name in the string (remoives file extensions)
			inputFileEnd=inputFileStart
			for char in inputName[inputFileStart:]:
				if char == ".":
					break
				else:
					inputFileEnd+=1


			# If the user has specified a file path, uses that instead of the default
			if not ammendedPath == "":

				outputFileName=ammendedPath
			else:

				outputFileName=(inputName[:inputFileStart])


			outputFileName +=  inputName[inputFileStart:inputFileEnd] + ammendedName + ".vcf.gz"
			outputNames.append(outputFileName)

	
	return outputNames

# Sorts arguments provided by the command line, and suplies defaults if no arguments are provided for some options. A complete set of parameters is returned as a dictionary
def SetParam(userInput):
	# Creates lists which contain the input out output file names
	inputSequences=[]
	# Creates two dictionaries, one which loads default parameters, while the other tracks if parameters were modified
	userParameters={"-d":8, "-D":95, "-f": "GRCh38_no_alt.fa", 
	"-l": False, "-o": "_filt", "input": inputSequences, "debuggingMode":False, "-p":"", "-q":5}
	modifiedParameters={"-d":False, "-D":False, "-f":False, "-l":False, "-o":False, "-p":False, "-q":False}

	# Cycles through each argument
	currentCom=-1
	optionAdded=False
	for args in userInput:

		currentCom +=1
		if optionAdded:
			optionAdded=False
			continue

		# If an input one of these option
		if (args  == "-f" or args == "-l" or args == "-d" or args == "-D" or args == "-o" or args == "-p"):

			# If there is is a paramter for the specified option
			if (currentCom+1) < len(userInput): 

				# Modifies Dictionary values if they input types match
				if ((args == "-f" or args == "-l" or args == "-o" or args == "-p") and isinstance(userInput[currentCom+1], str)):

					userParameters[args]=userInput[currentCom+1]
					modifiedParameters[args]=True
					optionAdded=True

				# 
				elif ((args == "-d" or args == "-D") and userInput[currentCom+1].isdigit()):
					
					userParameters[args]=userInput[currentCom+1]
					modifiedParameters[args]=True
					optionAdded=True

				else:
					#Replace this with an expection
					print("The specified paramer for " + args + " is invlid")
					exit()

			else:
				print("No parameter was specified for " + args)
				exit()

		# Activates debugging mode
		elif args == "-debug":

			userParameters["debuggingMode"]=True


		elif isinstance(args, str):

			inputSequences.append(args)
			modifiedParameters["seqInput"]=True

		else:

			print("Unknown input: " + args + ", ignoring...")


	while not os.path.isfile(userParameters["-f"]):

		print("The file " + userParameters["-f"] + " could not be found")
		exit()
		
	while not os.path.isfile(userParameters["-l"]):

		print("The file " + str(userParameters["-l"]) + " could not be found")
		exit()
		

	userParameters["output"]=setOutputFileNames(inputSequences, userParameters["-p"], userParameters["-o"])

	print("")

	# Informs the user which options are being left at default
	for key, modified in modifiedParameters.items():
		if not modified:

			if key == "input" or key == "-l":
				continue

			else:
				print("[" + key + " was not supplied, using " + str(userParameters[key]) + " ]")

	return userParameters 


def getPercentile(fileName, percentile, minDepth=0):


	def getDepth(infoColumn):

		startIndex = 0;
		endIndex=0;

		# Filters out "DP" and the other info in the information column
		for i in range(startIndex, len(infoColumn)):
		
			if infoColumn[i:i+2] == "DP":
				startIndex = i
				break

		for i in range(startIndex, len(infoColumn)):
		
			if infoColumn[i:i+1] == ";":
				endIndex = i
				break

		return int(infoColumn[startIndex+3:endIndex])



	print()
	print("Getting depth\n")
	vcfFile = open(fileName)
	fileLine = vcfFile.readline()
	fileElements = ["#None"]

	# Loop through the file until the VC header is reached
	while "#CHROM" not in fileElements[0]:

		fileLine = vcfFile.readline()
		fileElements = fileLine.split()

	#Skips down one line to the start of the variant calls
	fileLine = vcfFile.readline()
	fileElements = fileLine.split()
	print(fileElements)
	print("")

	variantDepth = []
	qualityAll = []

	# Loops through each varaint call
	for variantLines in vcfFile:

		# Makes a list of each each element in the line
		qualityAll.append(float(variantLines.split()[5]))
		infoColumn = variantLines.split()[7]

		depth = getDepth(infoColumn)

		# IF the read depth falls below the threshold, it is not included in list of read depth
		if depth < minDepth: continue

		# Adds this read depth to the list
		variantDepth.append(depth)

	# Prints out the designated percentile of the depth and quality
	return [numpy.percentile(variantDepth, float(percentile)), numpy.percentile(qualityAll, float(100-percentile))]

def runProcesses(parameters):


	mpileupArgs = ["samtools", "mpileup", "-ugf", parameters["-f"]]
	if parameters["-l"]:
		mpileupArgs.extend(["-l", parameters["-l"]])

	bcfCallArgs = ["bcftools", "call", "-vmO", "v"]

	counter=1
	for file in parameters["input"]:

		tempFile = "bufferFile" + str(counter) + ".vcf"

		bcfCallArgs.extend(["-o", tempFile, "-"])
		mpileupArgs.append(file)
		mpileup = subprocess.Popen(mpileupArgs, stdout=subprocess.PIPE)

		bcfCall = subprocess.Popen(bcfCallArgs, stdin=mpileup.stdout)
		bcfCall.communicate(mpileup.stdout)
		
		# Max depth is stored in limitingParam[0], Quality is stored in limitingParam[1]
		limitingParam = getPercentile(tempFile, parameters["-D"], parameters["-d"])
		print("Using depth " + str(limitingParam[0]) + " and quality of " + str(limitingParam[1]))

		filterArgs = "DP>" + str(parameters["-d"]) + " && DP<" + str(limitingParam[0]) + " && QUAL>" + str(limitingParam[1])
		bcfFilterArgs = ["bcftools", "filter", "-i", filterArgs, "-O", "v", "-o", parameters["output"][0], tempFile]
		bcfFilter = subprocess.run(bcfFilterArgs)
		counter += 1



cmdArgs=GetStartingVariables()

parameters=SetParam(cmdArgs)


# If debuggind mode is enabled, print out the options and filenames listed in the dictionary
if parameters["debuggingMode"]: printRawParameters(parameters)


print("\nSequencaRescuer: Starting mpileup\n")

runProcesses(parameters)
