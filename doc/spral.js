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
       scrollItems = menuItems.map(function(){
          var hr = $(this).attr("href")
          var item = $(jqescape(hr));
          if(item.length) { return item; }
          });
   $(window).scroll(function(){
      var fromTop = $(this).scrollTop();

      // Build list of items that are off top of page
      var cur = scrollItems.map(function() {
            if($(this).offset().top < fromTop) return this;
         });
      // Select last one
      cur = cur[cur.length-1];
      var id = cur && cur.length ? cur[0].id : "";

      if(lastId !== id) {
         lastId = id;
         var mp = menuItems.parent();
         var mp2 = mp.removeClass("active");
         var mp3 = mp2.end();
         var mp4 = mp3.filter("[href=#"+jqescape(id)+"]");
         var mp5 = mp4.parent();
         mp5.addClass("active");
      }
   });

}
