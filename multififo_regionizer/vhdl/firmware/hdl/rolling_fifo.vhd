library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library unisim;
use unisim.vcomponents.all;

use work.regionizer_data.all;

entity rolling_fifo is
    --generic(
    --    FIFO_INDEX : natural := 0
    --);
    port(
        ap_clk   : in std_logic;
        d_in     : in particle;
        write_in : in std_logic;
        d_out    : out particle;
        valid_out : out std_logic;
        full      : in std_logic;
    -- debug
        dbg_w64         : out std_logic_vector(63 downto 0);
    -- end debug
        roll      : in  std_logic;
        roll_out  : out std_logic
    );
end rolling_fifo;

architecture Behavioral of rolling_fifo is
    signal d64, q64 : std_logic_vector(63 downto 0);
    signal raddr, waddr : std_logic_vector(14 downto 0);
    signal rptr : unsigned(5 downto 0) := (0=>'1', others => '0'); -- need to count up to 63
    signal wptr : unsigned(5 downto 0) := (others => '0');         -- need to count up to 63
    signal wren, valid_next : std_logic := '0';
    signal roll_delay: std_logic_vector(2 downto 0) := (others => '0');
    signal cache : particle;
    signal cache_valid, mem_out_valid : std_logic := '0';
    signal use_cache   : std_logic := '0';
begin
    ram: RAMB36E2
        generic map(
            CLOCK_DOMAINS => "COMMON",
            DOA_REG => 0,
            DOB_REG => 0,
            READ_WIDTH_A => 72,
            WRITE_WIDTH_B => 72,
            WRITE_MODE_A => "READ_FIRST", --"WRITE_FIRST",
            WRITE_MODE_B => "READ_FIRST" --"WRITE_FIRST"
        )
        port map(
            ADDRENA => '1',
            ADDRENB => '1',
            ADDRARDADDR => raddr,
            ADDRBWRADDR => waddr,
            CLKARDCLK => ap_clk,
            CLKBWRCLK => ap_clk,
            DINADIN => d64(31 downto 0),
            DINBDIN => d64(63 downto 32),
            DINPADINP => (others => '0'),
            DINPBDINP => (others => '0'),
            DOUTADOUT => q64(31 downto 0),
            DOUTBDOUT => q64(63 downto 32),
            DOUTPADOUTP => open,
            DOUTPBDOUTP => open,
            ENARDEN => valid_next, -- avoid collisions 
            ENBWREN => wren,
            REGCEAREGCE => '1',
            REGCEB => '0',
            RSTRAMARSTRAM => '0',
            RSTRAMB => '0',
            RSTREGARSTREG => '0',
            RSTREGB => '0',
            WEA => "0000",
            WEBWE => "11111111",
            CASDIMUXA => '0',
            CASDIMUXB => '0',
            CASDOMUXA => '0',
            CASDOMUXB => '0',
            CASDINA => (others => '0'),
            CASDINB => (others => '0'),
            CASDINPA => (others => '0'),
            CASDINPB => (others => '0'),
            CASDOMUXEN_A => '1',
            CASDOMUXEN_B => '1',
            CASINSBITERR => '0',
            CASINDBITERR => '0',
            CASOREGIMUXEN_A => '1',
            CASOREGIMUXEN_B => '1',
            CASOREGIMUXA => '0',
            CASOREGIMUXB => '0',
            INJECTSBITERR => '0',
            INJECTDBITERR => '0',
            ECCPIPECE => '1',
            SLEEP => '0'

        );

     --- combinatorical
     raddr(14 downto 12) <= (others => '0');
     waddr(14 downto 12) <= (others => '0');
     raddr(11 downto 6) <= std_logic_vector(rptr) when valid_next = '1' else (others => '0'); -- avoid collision
     waddr(11 downto 6) <= std_logic_vector(wptr);
     raddr(5 downto 0) <= (others => '0');
     waddr(5 downto 0) <= (others => '0');

     out_switch: process(cache,q64,use_cache)
     begin
         if use_cache = '1' then
             d_out <= cache;
         else
             d_out <= w64_to_particle(q64);
         end if;
      end process;
     valid_out_switch: process(cache_valid,mem_out_valid,use_cache)
     begin
         if use_cache = '1' then
             valid_out  <= cache_valid;
         else
             valid_out  <= mem_out_valid;
         end if;
      end process;


     roll_out <= roll_delay(2);
     

     logic: process(ap_clk) 
           variable rptr_next : unsigned(5 downto 0);
           variable full_and_valid_out : boolean;
        begin
            if rising_edge(ap_clk) then
                
  -- assert rptr /= wptr or wren = '0'  report "Unexpected collision at " & 
  --        " FIFO index " & integer'image(FIFO_INDEX) &
  --        " rptr = wptr = " & integer'image(to_integer(rptr)) &
  --        " wren = " & std_logic'image(wren) & 
  --        " valid_next = " & std_logic'image(valid_next) & 
  --        " d_in(pt = " & integer'image(to_integer(d_in.pt)) &
  --        " eta = " & integer'image(to_integer(d_in.eta)) &
  --        " phi = " & integer'image(to_integer(d_in.phi)) &
  --        " rest = " & integer'image(to_integer(d_in.rest)) &
  --        ") write_in = " & std_logic'image(write_in) & 
  --        ", roll = " & std_logic'image(roll)  severity warning;

                if roll = '1' then
                    wptr <= (0 => write_in, others => '0');
                elsif write_in = '1' then
                    wptr <= wptr + to_unsigned(1, wptr'length);
                end if;
                wren <= write_in;

                d64 <= particle_to_w64(d_in);

                roll_delay(0) <= roll;
                roll_delay(2 downto 1) <= roll_delay(1 downto 0);

                full_and_valid_out := roll_delay(1) = '0' and roll_delay(2) = '0' and full = '1' and ((use_cache = '1' and cache_valid = '1') or (use_cache = '0' and mem_out_valid = '1'));

                if roll_delay(0) = '1' then
                    rptr_next := to_unsigned(1, wptr'length);
                elsif not(full_and_valid_out) and valid_next = '1' then
                    rptr_next := rptr + to_unsigned(1, rptr'length);
                else
                    rptr_next := rptr;
                end if;
                rptr <= rptr_next;

                if rptr_next <= wptr then
                    valid_next <= '1';
                else
                    valid_next <= '0';
                end if;

                if full_and_valid_out then
                    if use_cache = '0' then
                        cache <= w64_to_particle(q64);
                        cache_valid <= mem_out_valid;
                    end if;
                    use_cache <= '1';
                else 
                    cache <= w64_to_particle(q64);
                    cache_valid <= mem_out_valid;
                    use_cache   <= '0';
                end if;
                mem_out_valid <= valid_next;

            end if;

        end process;


       dbg_w64(15 downto 0) <= (0 => valid_next, 1 => use_cache, 2 => cache_valid, 4 => wren, others => '0');
       dbg_w64(21 downto 16) <= std_logic_vector(rptr);
       dbg_w64(31 downto 22) <= (others => '0');
       dbg_w64(37 downto 32) <= std_logic_vector(wptr);
       dbg_w64(47 downto 38) <= (others => '0');
       dbg_w64(63 downto 48) <= std_logic_vector(cache.pt);
    
end Behavioral;
