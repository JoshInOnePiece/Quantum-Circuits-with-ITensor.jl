using ITensors
using Statistics
using Printf
include(raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\Quantum Circuits Library\FILE IO Quantum Circuits Library\fileBitAdder.jl")
include(raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\Quantum Circuits Library\FILE IO Quantum Circuits Library\fileTEBDCircuit.jl")
include(raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\Quantum Circuits Library\FILE IO Quantum Circuits Library\fileVariationAnsatz.jl")
let
    #Names of input and output files
    inputAdderFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\bitAdderInput.txt"
    outputAdderFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\bitAdderOutput.txt"
    inputVariationalFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\variationalAnsantzInput.txt"
    outputVariationalFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\variationalAnsantzOutput.txt"
    inputTEBDFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\tebdInput.txt"
    outputTEBDFileName = raw"C:\Users\manik\Documents\Quantum-Circuits-with-ITensor.jl\IO FILES\tebdOutput.txt"

    #opens files from which input and output is stored
    inputAdderFile = open(inputAdderFileName,"r")
    outputAdderFile = open(outputAdderFileName, "w")
    inputVariationalFile = open(inputVariationalFileName,"r")
    outputVariationalFile = open(outputVariationalFileName, "w")
    inputTEBDFile = open(inputTEBDFileName,"r")
    outputTEBDFile = open(outputTEBDFileName, "w")

    #Calling the Binary Adder Circuit
    fileAdderCircuit(inputAdderFile, outputAdderFile)
    fileVariationalCircuit(inputVariationalFile, outputVariationalFile)
    fileTEBDCircuit(inputTEBDFile, outputTEBDFile)

    #Closes files
    close(inputAdderFile)
    close(inputVariationalFile)
    close(inputTEBDFile)
    close(outputAdderFile)
    close(outputVariationalFile)
    close(outputTEBDFile)
end