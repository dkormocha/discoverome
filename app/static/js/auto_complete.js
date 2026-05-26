document.addEventListener("DOMContentLoaded", () => {
  $("#auto").autocomplete({
    source: function(request, response) {
      $.ajax({
        url: "/autocomplete",
        data: {
          term: request.term
        },
        success: response
      });
    },
    minLength: 3
  });
  $("#target").autocomplete({
    source: function(request, response) {
      $.ajax({
        url: "/path_autocomplete",
        data: {
          path: request.term
        },
        success: response
      });
    },
    minLength: 3
  });

});