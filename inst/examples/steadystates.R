\dontrun{
reactions <- eqnlist()
reactions <- addReaction(reactions, "Tca_buffer", "Tca_cyto", 
                         "import_Tca*Tca_buffer", "Basolateral uptake")
reactions <- addReaction(reactions, "Tca_cyto", "Tca_buffer", 
                         "export_Tca_baso*Tca_cyto", "Basolateral efflux")
reactions <- addReaction(reactions, "Tca_cyto", "Tca_canalicular", 
                         "export_Tca_cana*Tca_cyto", "Canalicular efflux")
reactions <- addReaction(reactions, "Tca_canalicular", "Tca_buffer", 
                         "transport_Tca*Tca_canalicular", "Transport bile")

mysteadies <- steadyStates(reactions)
print(mysteadies)
}
