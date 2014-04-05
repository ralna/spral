function jqescape( arg ) {
   return arg.replace( /(:|\.|\[|\])/g, "\\$1" );
}
function spralDocReady() {
   $('.TOC ul ul').hide(); 
   $('.toc_chapter').click(function() { 
      $(this).find('ul').slideToggle(); 
      $(this).find('.toc_indicator').toggleClass('toc_dArrow'); 
      $(this).find('.toc_indicator').toggleClass('toc_rArrow'); 
   }); 

   var lastId,
       tocMenu = $("#tocMasterUL"),
       menuItems = tocMenu.find("a"),
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
