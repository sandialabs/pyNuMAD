class Mesh():
    def __init__(self, nodes, elements, sets, sections, adhesiveNds, adhesiveEls, adhesiveElSet):
        """_summary_

        Args:
            nodes (_type_): _description_
            elements (_type_): _description_
            sets (_type_): _description_
            sections (_type_): _description_
            adhesiveNds (_type_): _description_
            adhesiveEls (_type_): _description_
            adhesiveElSet (_type_): _description_
        """
        self.nodes = nodes
        self.elements = elements
        self.sets = sets
        self.sections = sections
        self.adhesiveNds = adhesiveNds
        self.adhesiveEls = adhesiveEls
        self.adhesiveElSet = adhesiveElSet
        