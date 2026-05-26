document.addEventListener("DOMContentLoaded", () => {

  function createSlider(id, options) {
    const el = document.getElementById(id);

    if (!el) return null;

    noUiSlider.create(el, {
      connect: true,
      tooltips: true,
      ...options
    });

    return el;
  }

  const slider = createSlider("slider", {
    start: [4, 20],
    range: { min: 0, max: 20 }
  });

  const sliderPocketVol = createSlider("slider-pocket_vol", {
    start: [350, 10000],
    range: { min: 350, max: 10000 }
  });

  const sliderDrugScore = createSlider("slider-drug_score", {
    start: [0, 1],
    range: { min: 0, max: 1 }
  });

  const sliderDistance = createSlider("slider-distance", {
    start: [0, 8],
    range: { min: 0, max: 8 }
  });

  const sliderAccess = createSlider("slider-access", {
    start: [0, 110],
    range: { min: 0, max: 110 }
  });

  const sliderPromis = createSlider("slider-promis", {
    start: [0, 6000],
    range: {
      min: 0,
      "50%": 100,
      "60%": 250,
      max: 6000
    },
    pips: {
      mode: "count",
      values: 3,
      density: 5
    }
  });

  window.toggleSlider = function(sliderId, checkboxId) {
    const slider = document.getElementById(sliderId);
    const checkbox = document.getElementById(checkboxId);

    slider.toggleAttribute("disabled", !checkbox.checked);
  };

  [
    ["slider", "crCheckbox"],
    ["slider-pocket_vol", "pocketVolCheckbox"],
    ["slider-drug_score", "drugScoreCheckbox"],
    ["slider-distance", "distanceCheckbox"],
    ["slider-access", "accessCheckbox"],
    ["slider-promis", "promisCheckbox"]
  ].forEach(([sliderId, checkboxId]) => {
    toggleSlider(sliderId, checkboxId);
  });

  $("#adv_search").on("click", () => {

    const params = new URLSearchParams();

    params.append(
      "res_type",
      $("#resTypeSelect").val()
    );

    function appendSliderValues(name, sliderEl) {

      if (!sliderEl || sliderEl.hasAttribute("disabled")) {
        params.append(`${name}_min`, "None");
        params.append(`${name}_max`, "None");
        return;
      }

      const [min, max] = sliderEl.noUiSlider.get();

      params.append(`${name}_min`, min);
      params.append(`${name}_max`, max);
    }

    appendSliderValues("cr", slider);
    appendSliderValues("pocket_vol", sliderPocketVol);
    appendSliderValues("drug_score", sliderDrugScore);
    appendSliderValues("distance", sliderDistance);
    appendSliderValues("access", sliderAccess);
    appendSliderValues("promis", sliderPromis);

    params.append(
      "proteinlist",
      $("#proteinlist").val().trim().replace(/\n/g, ",")
    );

    params.append(
      "target",
      $("#target").val()
    );

    window.location.href = `/advance_search?${params.toString()}`;
  });

});