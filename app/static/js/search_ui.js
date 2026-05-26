document.addEventListener("DOMContentLoaded", () => {

  const advBtn =
    document.getElementById("adv-search-btn");

  const advPanel =
    document.getElementById("adv-search");

  if (advBtn && advPanel) {

    advBtn.addEventListener("click", () => {

      const isOpen =
        advPanel.style.display === "block";

      advPanel.style.display =
        isOpen ? "none" : "block";

      advBtn.classList.toggle("open", !isOpen);
    });
  }

  window.setMode = function(mode) {

    const proteinForm =
      document.getElementById("proteinForm");

    const pathwayForm =
      document.getElementById("targetForm");

    const proteinBtn =
      document.getElementById("btn-protein");

    const pathwayBtn =
      document.getElementById("btn-pathway");

    if (mode === "protein") {

      proteinForm.style.display = "";
      pathwayForm.style.display = "none";

      proteinBtn.classList.add("active");
      pathwayBtn.classList.remove("active");

    } else {

      pathwayForm.style.display = "";
      proteinForm.style.display = "none";

      pathwayBtn.classList.add("active");
      proteinBtn.classList.remove("active");
    }
  };

});