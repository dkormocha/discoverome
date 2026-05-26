document.addEventListener("DOMContentLoaded", () => {

  if ($("#list").length) {

    new DataTable("#list", {

      pageLength: 10,

      language: {
        search: "",
        searchPlaceholder: "Filter pathways…",
        lengthMenu: "Show _MENU_ entries"
      }
    });
  }

});