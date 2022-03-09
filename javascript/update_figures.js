
function update_in_vitro_decays_figure() {
  var filename = "./figures/figS01_in_vitro_decays/2021-11-18_PBS_35C_decays_" + document.getElementById("raw_decay_pH_yaxis").value; 
  var image = document.getElementById("raw_decay_pH");
  image.src = filename;
  }  