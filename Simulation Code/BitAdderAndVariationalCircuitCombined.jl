using ITensors
using Statistics
using Printf
using Plots
include(raw"Quantum Circuits Library\FILE IO Quantum Circuits Library\fileBitAdder.jl")
include(raw"Quantum Circuits Library\FILE IO Quantum Circuits Library\fileTEBDCircuit.jl")
include(raw"Quantum Circuits Library\FILE IO Quantum Circuits Library\fileVariationAnsatz.jl")
let
    #Names of input and output files
    print("Imported Libraries\n")
    print("Opening Files\n")
    inputAdderFileName = raw"IO FILES\bitAdderInput.txt"
    outputAdderFileName = raw"IO FILES\bitAdderOutput.txt"
    inputVariationalFileName = raw"IO FILES\variationalAnsantzInput.txt"
    outputVariationalFileName = raw"IO FILES\variationalAnsantzOutput.txt"
    inputTEBDFileName = raw"IO FILES\tebdInput.txt"
    outputTEBDFileName = raw"IO FILES\tebdOutput.txt"

    #opens files from which input and output is stored
    inputAdderFile = open(inputAdderFileName,"r")
    outputAdderFile = open(outputAdderFileName, "w")
    inputVariationalFile = open(inputVariationalFileName,"r")
    outputVariationalFile = open(outputVariationalFileName, "w")
    inputTEBDFile = open(inputTEBDFileName,"r")
    outputTEBDFile = open(outputTEBDFileName, "w")

    print("Starting Circuit Simulations\n")
    #Calling the Binary Adder Circuit
    fileAdderCircuit(inputAdderFile, outputAdderFile)
    print("Finished Adder Circuit\n")
    fileVariationalCircuit(inputVariationalFile, outputVariationalFile)
    print("Finished Variational Circuit\n")
    fileTEBDCircuit(inputTEBDFile, outputTEBDFile)
    print("Finished TEBD Circuit\n")

    print("Closing files\n")
    #Closes files
    close(inputAdderFile)
    close(inputVariationalFile)
    close(inputTEBDFile)
    close(outputAdderFile)
    close(outputVariationalFile)
    close(outputTEBDFile)
    print("Simulation Complete")
end