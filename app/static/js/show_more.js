$(document).ready(function(){
    
    $('#more').click(function(){
        $('#text').removeClass("less");
        $("#less").show();
        $(this).hide();
    });
    
    
    $('#less').click(function(){
        $('#text').addClass("less");
        $(this).hide();
        $("#more").show();
    });
})