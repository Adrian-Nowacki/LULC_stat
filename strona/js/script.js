function menu_show() {
    var x = document.getElementById("navDemo");
    if (x.className.indexOf("w3-show") == -1) {
      x.className += " w3-show";
    } else { 
      x.className = x.className.replace(" w3-show", "");
    }
  }

  
  function replaceClass(id, oldClass, newClass) {
    var elem = $(`#${id}`);
    if (elem.hasClass(oldClass)) {
        elem.removeClass(oldClass);
    }
    elem.addClass(newClass);
}

  $(document).ready(function(){       
    var scroll_start = 0;
    var startchange = $('.w3-content');
    var offset = startchange.offset();
    $(document).scroll(function() { 
       scroll_start = $(this).scrollTop();
       if(scroll_start > offset.top) {
           $('.navbar').css({
    'font-size' : '0.90rem',
    'background-color': '#333333',
    'transition': 'font-size ease 1s'
 });

 replaceClass("menu_kafelkow", "w3-large", "");
        } else {
           $('.navbar').css({
             'background-color': 'transparent',
             'font-size': '1.2rem'
           });
        }
    });
 });
 
 
 
 let scroll_procent = () => {
     let scrollProgress = document.getElementById("procent");
     let progressValue = document.getElementById("procent-value");
     let pos = document.documentElement.scrollTop;
     let calcHeight = document.documentElement.scrollHeight - document.documentElement.clientHeight;
     let scrollValue = Math.round( pos * 100 / calcHeight);
     scrollProgress.style.background = `conic-gradient(#1e9477 ${scrollValue}%, #b6e2d7 ${scrollValue}%)`;
     progressValue.textContent = `${scrollValue}%`;
 }
 window.onscroll = scroll_procent;
 window.onload = scroll_procent;
 
 