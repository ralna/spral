function jqescape( arg ) {
   return arg.replace( /(:|\.|\[|\])/g, "\\$1" );
}
function spralDocReady() {
   // Hide all entries below "chapter" in TOC, but allow user to toggle their
   // expansion
   $('.TOC ul ul').hide(); 
   $('.toc_indicator').click(function() { 
      $(this).parent().find('ul').slideToggle(); 
      $(this).parent().find('.toc_indicator').toggleClass('toc_dArrow'); 
      $(this).parent().find('.toc_indicator').toggleClass('toc_rArrow'); 
   }); 

   // Setup menu variables
   var tocMenu = $("#tocMasterUL"),
       menuItems = tocMenu.find("a");

   // Unhide current "chapter"
   var chapterid = $('#pkgHeader h1 a').first().attr("id");
   var chul = menuItems.filter("[href=#"+jqescape(chapterid)+"]").closest("li");
   chul.find('.toc_indicator').addClass('toc_dArrow');
   chul.find('.toc_indicator').removeClass('toc_rArrow');
   chul.find("ul").show();

   // Highlight current section with .active class as user scrolls.
   var lastId,
       headerHeight = $("#pgHeader").outerHeight(),
       scrollItems = menuItems.map(function(){
          var hr = $(this).attr("href")
          var item = $(jqescape(hr));
          if(item.length) { return item; }
          });
   $(window).scroll(function(){
      var fromTop = $(this).scrollTop()+headerHeight;

      // Build list of items that are off top of page
      var cur = scrollItems.map(function() {
            if($(this).offset().top < fromTop) return this;
         });
      // Select last one
      cur = cur[cur.length-1];
      var id = cur && cur.length ? cur[0].id : "";

      if(lastId !== id) {
         lastId = id;
         menuItems.parent().removeClass("active")
                  .end().filter("[href=#"+jqescape(id)+"]").parent().addClass("active");
      }
   });

}
