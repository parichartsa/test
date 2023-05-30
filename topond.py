class PhylogeneticTreeNode:
    def __init__(self, sequence):
        self.sequence = sequence
        self.leftChild = None
        self.rightChild = None

class PhylogeneticTreeConstruction:
    @staticmethod
    def constructPhylogeneticTree(sequences, distanceMatrix):
        numSequences = len(sequences)

        # Create leaf nodes for each sequence
        leafNodes = []
        for sequence in sequences:
            leafNode = PhylogeneticTreeNode(sequence)
            leafNodes.append(leafNode)

        while len(leafNodes) > 1:
            # Find the pair of nodes with the smallest distance
            minDistance = float('inf')
            minI = 0
            minJ = 0
            for i in range(len(leafNodes)):
                for j in range(i + 1, len(leafNodes)):
                    distance = distanceMatrix[i][j]
                    if distance < minDistance:
                        minDistance = distance
                        minI = i
                        minJ = j

            # Create a new internal node and set the pair as its children
            newNode = PhylogeneticTreeNode("")
            newNode.leftChild = leafNodes[minI]
            newNode.rightChild = leafNodes[minJ]

            # Remove the pair from leafNodes and add the new node
            leafNodes.remove(leafNodes[minJ])
            leafNodes.remove(leafNodes[minI])
            leafNodes.append(newNode)

            # Update the distance matrix
            updatedDistanceMatrix = [[0] * (numSequences - 1) for _ in range(numSequences - 1)]
            newRow = 0
            for i in range(numSequences):
                if i != minI and i != minJ:
                    newCol = 0
                    for j in range(numSequences):
                        if j != minI and j != minJ:
                            updatedDistanceMatrix[newRow][newCol] = distanceMatrix[i][j]
                            newCol += 1
                    newRow += 1

            # Calculate new distances between the new node and remaining nodes
            for i in range(numSequences):
                if i != minI and i != minJ:
                    dist1 = distanceMatrix[minI][i]
                    dist2 = distanceMatrix[minJ][i]
                    newDistance = (dist1 + dist2 - minDistance) // 2
                    updatedDistanceMatrix[newRow][i if i < minI else i - 1] = newDistance
                    updatedDistanceMatrix[i if i < minI else i - 1][newRow] = newDistance

            distanceMatrix = updatedDistanceMatrix
            numSequences -= 1

        return leafNodes[0]  # Return the root of the constructed tree

    @staticmethod
    def printPhylogeneticTree(root, indent):
        if root:
            print(indent + "- " + root.sequence)
            PhylogeneticTreeConstruction.printPhylogeneticTree(root.leftChild, indent + "  |")
            PhylogeneticTreeConstruction.printPhylogeneticTree(root.rightChild, indent + "  |")

if __name__ == "__main__":
    sequences = [
        "T--TACTCCACACAC",
        "TGCTACTGCAGACAC",
        "TACTACTGACCCAC",
        "TACCACTCCACACCC"
    ]

    distanceMatrix = [
        [0, 11, 12, 11],
        [0, 0, 12, 13],
        [0, 0, 0, 11],
        [0, 0, 0, 0]
    ]

    root = PhylogeneticTreeConstruction.constructPhylogeneticTree(sequences, distanceMatrix)
    PhylogeneticTreeConstruction.printPhylogeneticTree(root, "")
